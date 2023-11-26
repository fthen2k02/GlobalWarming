#!/bin/bash

shopt -s nullglob

if ! test -d "./todo"; then
	mkdir todo
fi

if ! test -d "./done"; then
	mkdir done
fi

> "./done/list.txt"

while true; do
	for wav_file in todo/*.wav; do
		echo "Processing ${wav_file}..."

		fixed_name=${wav_file//：/-}
		fixed_name=${fixed_name//:/-}
		fixed_name=${fixed_name//｜/-}
		fixed_name=${fixed_name//|/-}
		fixed_name=${fixed_name//＂/-}
		fixed_name=${fixed_name//\"/-}
		fixed_name=${fixed_name//⧸/-}

		if test "${wav_file}" != "${fixed_name}"; then
			echo "Invalid file name: \"${wav_file}\". Renaming to \"${fixed_name}\"..."
			mv "${wav_file}" "${fixed_name}"
		fi

		rm -f solution_*.txt
		if test -f results.txt; then
			cp results.txt results_bak.txt
		fi

		if ! ./process_file.exe "${fixed_name}" results.txt -1 10 2000 0; then
			echo "File processing returned an error. Aborting."
			exit 1
		fi
		#mv "${fixed_name}" "./done"
		rm "${fixed_name}"
		echo "${fixed_name}" >> "./done/list.txt"
		sync
	done

	sleep 2
done