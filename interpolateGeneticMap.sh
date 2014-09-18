# usage: THIS [map file]

cat $1 | awk -v i=$1 '
BEGIN{
	count=0;
	while(getline < i) {
		if ($3!=0) {
			if (count++ > 0) {
				map[count] = ($3-lastCm)/($4-lastMb);
				lastCm = $3;
				lastMb = $4;
			}
			else {
                                firstCm = $3;
				firstMb = $4;
				lastCm = $3;
				lastMb = $4;
			}
		}
	}
	map[1]=map[2];
	map[count+1]=map[count];
	count = 1;
}

$3!=0{
      	count++;
        printf("%s\t%s\t%.11f\t%s\n", $1, $2, $3, $4);
        lastMb = $4;
        lastCm = $3;
}

$3==0 {

# if not passed first yet (count == 1) compute distance to first and use that + rate to get map
	rate = map[count];
	if (count == 1) {
		dist = firstCm-(firstMb-$4)*rate;
		$3 = dist;
		printf("%s\t%s\t%.10f\t%s\n", $1, $2, $3, $4);
		lastMb=$4;
		lastCm=$3;
	}
# else use last rate
	else {
		dist = lastCm + rate*($4-lastMb);
		$3 = dist;
		printf("%s\t%s\t%.10f\t%s\n", $1, $2, $3, $4);
                lastMb=$4;
                lastCm=$3;
	}
}'
