#!/usr/local/bin/bash

#script to create a html file for each sample with qcgrind plots and stats

cat << EOF
<html><body><h1>Bamcheck output for: $1</h1>
<hr>
<h2> acgt-cycles </h2>
<img src=$2 >
<hr>
<h2> coverage </h2>
<img src=$3 >
<hr>
<h2> gc-content </h2>
<img src=$4 >
<hr>
<h2> gc-depth </h2>
<img src=$5 >
<hr>
<h2> indel-cycles </h2>
<img src=$6 >
<hr>
<h2> indel-dist </h2>
<img src=$7 >
<hr>
<h2> insert-size </h2>
<img src=$8 >
<hr>
<h2> quals </h2>
<img src=$9 >
<hr>
<h2> quals2 </h2>
<img src=${10} >
<hr>
<h2> quals3 </h2>
<img src=${11} >
<hr>
<h2> quals-hm </h2>
<img src=${12} >

</html></body>
EOF