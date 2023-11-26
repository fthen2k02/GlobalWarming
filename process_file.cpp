/*
 * MIT License (https://opensource.org/license/mit/)
 *
 * Copyright (C) 2023 The Infinity Chiller
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <limits>
#include <map>
#include <chrono>
#include <iomanip>
#include <string>
#include <cctype>
#include <list>
#include <mutex>
#include <condition_variable>

#include "AudioFile.h"

#define QUASIPERIOD_MIN         0.0004    // s
#define QUASIPERIOD_MAX         0.014     // s

/*
 * To reduce the variation of quasiperiods in the same segmentation, the range of allowed quasiperiods will be divided
 * into this number of subintervals which will be tested independently and whose limits will form a geometric
 * progression. Any interval [some_quasiperiod, maxQuasiperiodRatio * some_quasiperiod] that is included in
 * [QUASIPERIOD_MIN, QUASIPERIOD_MAX] will be included in exactly one of these subintervals.
 * A grater value for NUM_SUBINTERVALS will reduce the risk of errors, but will impact the execution time linearly
 * unless parallelized.
 */
#define NUM_SUBINTERVALS        20

using namespace std;

struct solData
{
    int subintervalNum, channel;
    size_t firstSegment;
    bool omitFirstMin;

    solData(int subintervalNum, int channel, size_t firstSegment, bool omitFirstMin) :
        subintervalNum(subintervalNum), channel(channel), firstSegment(firstSegment), omitFirstMin(omitFirstMin) { }
};

struct threadData
{
    thread t;
    vector<int> segmentation;
    multimap<double, solData> solutions;
    pair<int, int> progress;
};

struct resultsFileEntry
{
    double error;
    int channel, startSample, length;
    string wavPath;
};

AudioFile<double> haystack;
int numChannels, numSamples;
double maxError;
vector<pair<int, int>> needleExtremes;
vector<pair<bool, bool>> needleExtremesKnown;
double maxQuasiperiodRatio; /* Maximum ratio between two "quasiperiods" within the searched-for waveform. Must be > 1.0. */
vector<vector<threadData>> threads;
int needleSum, maxSolsFromInput, maxSolsInOutput;
bool maxErrorDefined, maxSolsFromInputDefined, maxSolsInOutputDefined;
bool createSolFiles;
long long fullProgress;
condition_variable cv;

void readNeedle()
{
    ifstream f("needle.txt");
    int max, min;

    f >> maxQuasiperiodRatio;

    needleSum = 0;
    while (f >> max >> min)
    {
        if (max > -1)
            needleSum += max;
        if (min > -1)
            needleSum += min;

        needleExtremes.push_back(make_pair(max, -min));
        needleExtremesKnown.push_back(make_pair(max > -1, min > -1));
    }
}

inline double compareAgainstPrevious(int channel, int pos, int length)
{
    double err = 0;

    for (int i = pos; i < pos + length; ++i)
    {
        auto sample1 = haystack.samples[channel][i];
        auto sample2 = haystack.samples[channel][i - length];

        auto diff = sample1 - sample2;
        auto coef = (abs(sample1) + abs(sample2)) / 2 * 3;
        coef = coef * coef + 1;

        err += diff * diff * coef;
    }

    return err / length;
}

inline bool eval(const list<pair<double, double>>& extremes, int firstMaxPos, int firstMinPos, double &error)
{
    double sum = 0;
    auto N = needleExtremes.size();
    int pos;

    pos = 0;
    for (const auto &u : extremes)
    {
        if (pos >= firstMaxPos && pos < firstMaxPos + N && needleExtremesKnown[pos - firstMaxPos].first)
            sum += abs(u.first);
        if (pos >= firstMinPos && pos < firstMinPos + N && needleExtremesKnown[pos - firstMinPos].second)
            sum += abs(u.second);
        ++pos;
    }

    if (sum < numeric_limits<float>::min())
        return false;

    double scale = needleSum / sum;

    error = 0;
    pos = 0;
    for (const auto &u : extremes)
    {
        if (pos >= firstMaxPos && pos < firstMaxPos + N && needleExtremesKnown[pos - firstMaxPos].first)
        {
            auto diff = u.first * scale - needleExtremes[pos - firstMaxPos].first;
            error += diff * diff;
        }
        if (pos >= firstMinPos && pos < firstMinPos + N && needleExtremesKnown[pos - firstMinPos].second)
        {
            auto diff = u.second * scale - needleExtremes[pos - firstMinPos].second;
            error += diff * diff;
        }
        ++pos;
    }

    return true;
}

void segment(int subintervalNum, int channel, int minQuasiPeriodSamples, int maxQuasiPeriodSamples)
{
    minQuasiPeriodSamples = max(minQuasiPeriodSamples, 1);
    maxQuasiPeriodSamples = max(maxQuasiPeriodSamples, 1);

    auto &myData = threads[subintervalNum][channel];
    auto &pos = myData.progress.first;

    if (numSamples >= minQuasiPeriodSamples * 2)
    {
        auto minErr = numeric_limits<double>::max();
        int bestQp;
        for (int qp = minQuasiPeriodSamples; qp <= maxQuasiPeriodSamples && 2 * qp <= numSamples; ++qp)
        {
            auto tempErr = compareAgainstPrevious(channel, qp, qp);
            if (tempErr < minErr)
            {
                minErr = tempErr;
                bestQp = qp;
            }
        }

        myData.segmentation.push_back(bestQp);
        pos = bestQp;

        while (pos + minQuasiPeriodSamples <= numSamples)
        {
            auto minErr = numeric_limits<double>::max();
            for (int qp = minQuasiPeriodSamples; qp <= maxQuasiPeriodSamples && pos + qp <= numSamples && qp <= pos; ++qp)
            {
                auto tempErr = compareAgainstPrevious(channel, pos, qp);
                if (tempErr < minErr)
                {
                    minErr = tempErr;
                    bestQp = qp;
                }
            }

            myData.segmentation.push_back(pos + bestQp);
            pos += bestQp;
        }
    }

    if (pos < numSamples)
    {
        myData.segmentation.push_back(numSamples);
        pos = numSamples;
    }

    list<pair<double, double>> extremes;
    for (unsigned int i = 0; i < myData.segmentation.size(); ++i)
    {
        int segmentStart = !i ? 0 : myData.segmentation[i - 1];
        double max, min;

        max = min = haystack.samples[channel][segmentStart];
        for (int j = segmentStart + 1; j < myData.segmentation[i]; ++j)
        {
            auto temp = haystack.samples[channel][j];
            if (temp > max)
                max = temp;
            if (temp < min)
                min = temp;
        }

        extremes.push_back(make_pair(max, min));

        if (extremes.size() == needleExtremes.size() + 1)
        {
            double tempError;
            if (eval(extremes, 0, 1, tempError) && (!maxErrorDefined || tempError < maxError))
                myData.solutions.insert(make_pair(tempError, solData(subintervalNum, channel, i - needleExtremes.size(), true)));
            extremes.pop_front();
        }

        if (extremes.size() == needleExtremes.size())
        {
            double tempError;
            if (eval(extremes, 0, 0, tempError) && (!maxErrorDefined || tempError < maxError))
                myData.solutions.insert(make_pair(tempError, solData(subintervalNum, channel, i - needleExtremes.size() + 1, false)));
        }

        if (maxSolsFromInputDefined)
            while (myData.solutions.size() > maxSolsFromInput)
                myData.solutions.erase(prev(myData.solutions.end()));

        myData.progress.second = myData.segmentation[i];
    }
}

void solInfoToFile(unsigned int solNum, solData& info)
{
    char fileName[30];
    sprintf_s(fileName, "solution_%d.txt", solNum);
    ofstream g(fileName);

    int add = info.omitFirstMin ? 1 : 0;

    for (unsigned int i = 0; i < needleExtremes.size() + add; ++i)
    {
        int start = info.firstSegment + i == 0 ? 0 : threads[info.subintervalNum][info.channel].segmentation[info.firstSegment + i - 1];
        int end = threads[info.subintervalNum][info.channel].segmentation[info.firstSegment + i];
        for (int j = start; j < end; ++j)
            g << (i % 2) << ' ' << haystack.samples[info.channel][j] << '\n';
    }
}

bool isDone()
{
    long long sum_segm = 0, sum_eval = 0;

    for (int i = 0; i < NUM_SUBINTERVALS; ++i)
        for (int j = 0; j < numChannels; ++j)
        {
            sum_segm += threads[i][j].progress.first;
            sum_eval += threads[i][j].progress.second;
        }

    cerr << '\r' << sum_segm * 100 / fullProgress << "% " << sum_eval * 100 / fullProgress << "%";
    return sum_segm == fullProgress && sum_eval == fullProgress;
}

void showProgress()
{
    fullProgress = 1LL * numSamples * numChannels * NUM_SUBINTERVALS;

    if (fullProgress)
    {
        mutex mtx;
        unique_lock<mutex> lck(mtx);
        while (!cv.wait_for(lck, chrono::milliseconds(500), isDone));
    }

    cerr << "\r100% 100%\n";
}

istream& operator>>(istream& f, resultsFileEntry& entry)
{
    f >> entry.error >> entry.channel >> entry.startSample >> entry.length;
    getline(f, entry.wavPath);

    unsigned int i = 0;
    while (i < entry.wavPath.size() && isspace(entry.wavPath[i]))
        ++i;
    entry.wavPath = entry.wavPath.substr(i);

    return f;
}

int main(int argc, char* argv[])
{
    if (argc != 7)
    {
        cerr << "Usage: " << argv[0] << " <wav_file> <results_file> <max_error> <max_solutions_from_input> <max_solutions_in_output> <create_sol_files>\n";
        return 1;
    }

    if (!haystack.load(argv[1]))
    {
        cerr << "Unable to load WAV data from " << argv[1] << ".\n";
        return 1;
    }

    int sampleRate = haystack.getSampleRate();
    numSamples = haystack.getNumSamplesPerChannel();
    numChannels = haystack.getNumChannels();

    sscanf_s(argv[3], "%lf", &maxError);
    maxErrorDefined = maxError >= 0;
    maxSolsFromInput = atoi(argv[4]);
    maxSolsFromInputDefined = maxSolsFromInput >= 0;
    maxSolsInOutput = atoi(argv[5]);
    maxSolsInOutputDefined = maxSolsInOutput >= 0;
    createSolFiles = atoi(argv[6]) != 0;

    readNeedle();

    threads.resize(NUM_SUBINTERVALS);
    for (int i = 0; i < NUM_SUBINTERVALS; ++i)
        threads[i].resize(numChannels);

    thread progressLogger(showProgress);

    /* The ratio between the limits of each subinterval will be maxQuasiperiodRatio + delta. */
    double delta = pow(pow(maxQuasiperiodRatio, NUM_SUBINTERVALS - 1) * QUASIPERIOD_MAX / QUASIPERIOD_MIN, 1.0l / NUM_SUBINTERVALS) - maxQuasiperiodRatio;
    double subintervalRatio = maxQuasiperiodRatio + delta;
    double lower = QUASIPERIOD_MIN;
    double upper = lower * subintervalRatio;
    int lowerSamples = lround(lower * sampleRate);
    int upperSamples = lround(upper * sampleRate);

    for (int i = 0; i < NUM_SUBINTERVALS; ++i)
    {
        for (int j = 0; j < numChannels; ++j)
            threads[i][j].t = thread(segment, i, j, lowerSamples, upperSamples);

        lower = upper / maxQuasiperiodRatio;
        upper = lower * subintervalRatio;
        lowerSamples = static_cast<int>(upperSamples / maxQuasiperiodRatio) + 1;
        upperSamples = lround(upper * sampleRate);
    }

    for (unsigned int i = 0; i < NUM_SUBINTERVALS; ++i)
        for (int j = 0; j < numChannels; ++j)
            threads[i][j].t.join();

    cv.notify_one();
    progressLogger.join();

    multimap<double, solData> all;
    for (unsigned int i = 0; i < NUM_SUBINTERVALS; ++i)
        for (int j = 0; j < numChannels; ++j)
        {
            all.merge(threads[i][j].solutions);
            threads[i][j].solutions.clear();
        }

    int i = 0;
    for (auto &u : all)
    {
        if (maxSolsFromInputDefined && i >= maxSolsFromInput)
            break;

        cout << i << ' ' << u.first << ' ' << u.second.channel << ' ' << u.second.omitFirstMin << ' ';
        cout << (u.second.firstSegment == 0 ? 0 : threads[u.second.subintervalNum][u.second.channel].segmentation[u.second.firstSegment - 1]) << '\n';

        if (createSolFiles)
            solInfoToFile(i, u.second);

        ++i;
    }

    ifstream resultsIn(argv[2]);
    vector<resultsFileEntry> results;
    if (resultsIn.good())
    {
        resultsFileEntry entry;

        while (resultsIn >> entry)
            results.push_back(entry);

        resultsIn.close();
    }

    ofstream resultsOut(argv[2]);
    auto u = all.begin();
    int iOldResults = 0, iNewResults = 0;
    bool hasOld = iOldResults < results.size();
    bool hasNew = u != all.end() && (!maxSolsFromInputDefined || iNewResults < maxSolsFromInput);

    while ((hasOld || hasNew) && (!maxSolsInOutputDefined || iOldResults + iNewResults < maxSolsInOutput))
    {
        bool chooseOld = hasOld && hasNew ? results[iOldResults].error < u->first : hasOld;

        if (chooseOld)
        {
            resultsOut << fixed << setprecision(10) << results[iOldResults].error << '\t' << results[iOldResults].channel << '\t';
            resultsOut << results[iOldResults].startSample << '\t' << results[iOldResults].length << '\t';
            resultsOut << results[iOldResults].wavPath << '\n';

            ++iOldResults;
            hasOld = iOldResults < results.size();
        }
        else
        {
            resultsOut << fixed << setprecision(10) << u->first << '\t' << u->second.channel << '\t';
            int startSample = u->second.firstSegment == 0 ? 0 : threads[u->second.subintervalNum][u->second.channel].segmentation[u->second.firstSegment - 1];
            resultsOut << startSample << '\t';
            auto lastSegment = u->second.firstSegment + needleExtremes.size();
            if (!u->second.omitFirstMin)
                --lastSegment;
            resultsOut << threads[u->second.subintervalNum][u->second.channel].segmentation[lastSegment] - startSample << '\t';
            resultsOut << argv[1] << '\n';

            ++u;
            ++iNewResults;
            hasNew = u != all.end() && (!maxSolsFromInputDefined || iNewResults < maxSolsFromInput);
        }
    }

    return 0;
}