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
	cores="$3"
	offset="\"$4\""
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
				cores="`echo "$data"|cut -f3|uniq|tac`"
				case "$benchmarkName" in
					"\"add 1 w/ large strides"*)
						echo 'set ylabel "Add Operations [GFlop/s]"'
						;;
					*)
						echo 'set ylabel "Aggregate Bandwidth [GB/s]"'
						;;
				esac
				for o in $offsets; do
					inlineData=""
					first=true
					echo "set title \"${benchmarkName//\"/}, Offset $o\""
					for c in $cores; do
						test `coreCount "$c"` -eq 1 && o=0 # also plot single-threaded results on pages with offsets
						subset="`sliceSubset "$data" "$benchmarkName" "$c" "$o"`"
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
#plot
#'<grep "bzero (1 GiB)' \
#   using 1:($4+$3):4 with filledcurve title "User CPU Load" axes x1y1, \
#'' using 1:($4) with filledcurve x1 title "System CPU Load" axes x1y1, \
#'' using 1:($6+$4+$3):($4+$3) with filledcurve title "Waiting CPU" axes x1y1, \
#'' using 1:(min(3995,(max($27, $28) - 2000))*0.025) with lines lw 3 lc rgb "#909090" title "Fan" axes x1y1, \
#'' using 1:(max($7, (int($23) % 129000))*0.001) with lines lw 3 lc rgb "#00A000" title "GPU" axes x1y2, \
#'' using 1:(max($25, $26)*0.001) with lines lw 3 lc rgb "#FF9040" title "coretemp" axes x1y2, \
#'' using 1:($20*0.001) with lines lw 3 lc rgb "#008080" title columnheader axes x1y2, \
#'' using 1:($13*0.001) with lines lw 3 lc rgb "#800000" title columnheader axes x1y2, \
#'' using 1:($14*0.001) with lines lw 3 lc rgb "#808000" title columnheader axes x1y2

exit 0

maxy=7
test -n "$2" && test $2 -gt 0 && maxy=$2

xlabel="Memory Position [GB]"

dir=`dirname $0`

testnames="`cat "$1"|awk '{ if(FNR>2 && index($1, "benchmark.name") != 2) print }'|cut -d: -f1|
cut -c2-|awk '{
if (! ($0 in seen)) {
	seen[$0] = 1
	print
}
}'`"

psfile="$1.ps"
echo -n > "$psfile"

awkscript=`mktemp`
trap "rm $awkscript" EXIT INT HUP TERM

cat > $awkscript <<EOF
BEGIN {
   getline version
   getline header
   outHeader[1] = "Benchmark"
   data[1, 1] = ""
   line = 0
   col = 1
   row = 1
   lastcpu = -1
   while(0 != getline) {
      split(\$0, a, /"?\\t"?/)
      name  = substr(a[1], 2)
      if(length(filter) == 0 || index(name, filter) == 1) {
         cpu   = a[3]
         value = a[valueIndex]

         if(cpu != lastcpu) {
            row = 1
            lastcpu = cpu
            ++col
            outHeader[col] = "CPU" cpu
         }
         data[row, 1] = name
         if(invert) {
            data[row, col] = 1.0 / value
         } else {
            data[row, col] = value
         }
         if(length(factor) > 0) {
            data[row, col] *= factor
         }
         ++row
      }
   }
   for(i = 1; i <= length(outHeader); ++i) {
      printf("\\"%s\\"\\t", outHeader[i])
   }
   print ""
   for(r = 1; r < row; ++r) {
      printf("\\"%s\\"\\t", data[r, 1]);
      for(c = 2; c <= col; ++c) {
         printf("%f\\t", data[r, c]);
      }
      print ""
   }
}
EOF

oldIFS="$IFS"
IFS="${IFS# }"
for title in ${testnames}; do
	lmaxy="[0:${maxy}e9]"
	IFS="$oldIFS"
	csv=`mktemp`
	echo -n "#" > "$csv"
	case "$title" in
	*latency*)
		gawk -f "$awkscript" "-vfilter=$title:" "-vvalueIndex=10" "-vinvert=1" "$1" >> "$csv"
		tmp=`cut -f2 "$csv"|sort -un|tail -n1|cut -d. -f1`
		tmp=$(($tmp/10*10+10))
		lmaxy="[0:$tmp]"
		ylabel="Latency [cycles/read]"
		;;
	*)
		gawk -f "$awkscript" "-vfilter=$title:" "-vvalueIndex=8" "-vinvert=0" "-vfactor=1e-9" "$1" >> "$csv"
		lmaxy="[0:]"
		ylabel="Throughput [GB/s]"
		;;
	esac
	head=`cat "$csv"|head -n1|cut -f2-`
	width=`cat "$csv"|wc -l`
	width=$(($width-2))
	commands=""
	pos=2
	for h in $head; do
		commands="${commands} \
			\"$csv\" using $pos with linespoints title $h,"
		pos=$(($pos+1))
	done
	commands="${commands%,}"
	gnuplot <<EOF
set grid y
set terminal postscript color
set output "$csv.ps"
set title "$1         $title"
set ylabel "$ylabel"
set xlabel "$xlabel"
set xrange [0:$width]
set yrange ${lmaxy}
plot $commands
EOF
	cat "$csv.ps" >> "$psfile"
	rm "$csv" "$csv.ps"
done
IFS="$oldIFS"

# The Postscript file contains lines like this now:
# %%Title: /tmp/tmp.Rtu83pEYfL.ps
# /Title (/tmp/tmp.Rtu83pEYfL.ps)
#
# We'd like to have nicer title
sed -i \
  -e 's/%%Title: .*$/%%Title: '"$1"'/' \
  -e 's/\/Title \(.*\)$/\/Title ('"$1"')/' \
  "$psfile"

if which ps2pdf14 > /dev/null 2>&1; then
	ps2pdf14 "$psfile" "$1.pdf"
	rm "$psfile"
	echo "-> $1.pdf"
else
	echo "Install ps2pdf14 to get PDF output"
	echo "-> $psfile"
fi

# vim: noet sw=4 ts=4
