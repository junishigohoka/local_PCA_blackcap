BEGIN{
	OFS="\t"
}
{
	l=$2-$1
	r=$3
	L+=100*l*r
	print $1,r*1e+8,L
}
