#!/bin/bash
function miseq2iseq () {
	# From 0 to 17 ==> 11 (,)
	sed -i -e '4~4s/\!/\,/g' $1
	sed -i -e '4~4s/\"/\,/g' $1
	sed -i -e '4~4s/\#/\,/g' $1
	sed -i -e '4~4s/\$/\,/g' $1
	sed -i -e '4~4s/\%/\,/g' $1
	sed -i -e '4~4s/\&/\,/g' $1
	sed -i -e "4~4s/'/\,/g" $1
	sed -i -e '4~4s/(/\,/g' $1
	sed -i -e '4~4s/)/\,/g' $1
	sed -i -e '4~4s/\*/\,/g' $1
	sed -i -e '4~4s/\+/\,/g' $1
	#sed -i -e '4~4s/\,/\,/g' $1
	sed -i -e '4~4s/\-/\,/g' $1
	sed -i -e '4~4s/\./\,/g' $1
	sed -i -e '4~4s/\//\,/g' $1
	sed -i -e '4~4s/0/\,/g' $1
	sed -i -e '4~4s/1/\,/g' $1
	sed -i -e '4~4s/2/\,/g' $1

	# From 18 to 29 ==> 25 (:)
	sed -i -e '4~4s/3/\:/g' $1
	sed -i -e '4~4s/4/\:/g' $1
	sed -i -e '4~4s/5/\:/g' $1
	sed -i -e '4~4s/6/\:/g' $1
	sed -i -e '4~4s/7/\:/g' $1
	sed -i -e '4~4s/8/\:/g' $1
	sed -i -e '4~4s/9/\:/g' $1
	#sed -i -e '4~4s/\:/\:/g' $1
	sed -i -e '4~4s/\;/\:/g' $1
	sed -i -e '4~4s/</\:/g' $1
	sed -i -e '4~4s/\=/\:/g' $1
	sed -i -e '4~4s/>/\:/g' $1

	# From 30 to 40 ==> 37 (F)
	sed -i -e '4~4s/\?/F/g' $1
	sed -i -e '4~4s/\@/F/g' $1
	sed -i -e '4~4s/A/F/g' $1
	sed -i -e '4~4s/B/F/g' $1
	sed -i -e '4~4s/C/F/g' $1
	sed -i -e '4~4s/D/F/g' $1
	sed -i -e '4~4s/E/F/g' $1
	#sed -i -e '4~4s/F/F/g' $1
	sed -i -e '4~4s/G/F/g' $1
	sed -i -e '4~4s/H/F/g' $1
	sed -i -e '4~4s/I/F/g' $1
}
