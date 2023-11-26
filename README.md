The *San Andreas Mercenaries* update for Grand Theft Auto Online introduced [four new mysterious symbols](https://imagizer.imageshack.com/img923/9506/8meqGh.png), all of them found both on the ["??? Tee"](https://static.wikia.nocookie.net/gtawiki/images/9/97/SpecialClothing-GTAOee-%3F%3F%3FTee.jpg/revision/latest/) and in a secret facility situated under Fort Zancudo, each symbol associated to a room. One of these symbols contains a waveform behind a globe.

The purpose of this tool is to identify that waveform (the "needle") in a given WAV file (the "haystack"), based on the assumption that it is a sound wave of a spoken message, as justified [here](https://www.reddit.com/r/chiliadmystery/comments/17lngkh/comment/k7gmspw/?utm_source=reddit&utm_medium=web2x&context=3). It doesn't use AI, but an ad-hoc algorithm that only takes into account the maximum peaks in each "quasiperiod", both in the positive and negative directions, as they are the easiest to detect.

The complexity is Θ(N x Δ + N / Q x nq) ≈ Θ(N x Δ), where:
 - N is the number of samples in the haystack;
 - Δ is the size of the range within which the quasiperiod can be (depends on the sample rate);
 - Q is the number of samples in the average expected quasiperiod;
 - nq is the number of quasiperiods of the needle (77 in our case).

To read the WAV file, I used the [AudioFile](https://github.com/adamstark/AudioFile) library by Adam Stark, which is licensed under the MIT license. See the `AudioFile` subfolder for more information.