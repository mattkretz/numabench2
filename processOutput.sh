#!/bin/bash
fatal() {
	echo "$1" 1>&2
	exit 1
}
debug() {
	echo "$1" 1>&2
}

if ! test -f "$1"; then
	fatal "You need to specify the datafile!" >&2
fi

if ! which gawk >/dev/null 2>&1; then
	fatal "you need to have 'gawk' in your PATH"
fi

if ! which gnuplot >/dev/null 2>&1; then
	fatal "you need to have 'gnuplot' in your PATH"
fi

input="$1"
output="${input%.*}.pdf"

sliceOneBenchmark() {
	cut -f1,3,4,5,10 "$input" | grep -G "$1" | cut -f2-
}

sliceSubset() {
	data="$1"
	cores="$2"
	offset="\"$3\""
	echo "$data" | grep -G '".*"'".${offset}.${cores}" | cut -f1,4 | sed 's/"//g'
}

coreCount() {
	echo ",$1"|tr -cd ,|wc -c
}

plots() {
	cut -f1 "$input" | sort -u | while read benchmarkName; do
		case "$benchmarkName" in
			Version*)
				[[ "$benchmarkName" != "Version 3" ]] && fatal "Only Version 3 data files can be processed. (${benchmarkName})"
				;;
			"\"benchmark.name\"")
				;;
			*)
				data="`sliceOneBenchmark "$benchmarkName"`"
				offsets="`echo "$data"|cut -f2|tr -d \\"|sort -nu`"
				cores="`echo "$data"|cut -f3|uniq`"
				singleCores=(`echo "$cores"|grep -v ,`)
				case "$benchmarkName" in
					"\"add 1 w/ large strides"*)
						echo 'set ylabel "Add Operations [GFlop/s]"'
						;;
					*)
						echo 'set ylabel "Aggregate Bandwidth [GB/s]"'
						;;
				esac
				for firstCore in "${singleCores[@]}"; do
					firstCore=${firstCore%\"}
					coresSubset=(`echo "$cores"|grep -E "^${firstCore}[,\"]"|tac`)
					for o in $offsets; do
						inlineData=""
						first=true
						echo "set title \"${benchmarkName//\"/}, Offset $o\""
						for c in "${coresSubset[@]}"; do
							test `coreCount "$c"` -eq 1 && o=0 # also plot single-threaded results on pages with offsets
							subset="`sliceSubset "$data" "$c" "$o"`"
							test -z "$subset" && continue
							inlineData="${inlineData}${subset}
e
"
							if $first; then
								first=false
								echo -n "plot "
							else
								echo -n ", "
							fi
							echo -n "'-' using 1:(\$2*1e-9) with linespoints title ${c}"
						done
						echo -e "\n${inlineData}"
					done
				done
				;;
		esac
	done
}

#plots; exit

gnuplot <<EOF
set style line  1 lc rgbcolor "#CCCCCC"
set grid y ls 1
set grid x ls 1
set style line  1 lc rgbcolor "#AF3737" lw 2
set style line  2 lc rgbcolor "#3737AF" lw 2
set style line  3 lc rgbcolor "#37AF37" lw 2
set style line  4 lc rgbcolor "#AFAF37" lw 2
set style line  5 lc rgbcolor "#AF37AF" lw 2
set style line  6 lc rgbcolor "#37AFAF" lw 2
set style increment user
#set ytics nomirror 10
set yrange [0 : *]
#set y2label "Temperatur [°C]"
#set y2tics nomirror 10
#set y2range [10 : 110]
set xlabel "Memory Location [N × GiB]"
set terminal pdf noenhanced font "Droid Sans,10" size 29.7cm,21cm
set output "$output"
set style fill solid 0.3 noborder
#set key top left Left reverse
#set samples 400
max(a, b) = a > b ? a : b
min(a, b) = a < b ? a : b
$(plots)
EOF

# vim: noet sw=4 ts=4
