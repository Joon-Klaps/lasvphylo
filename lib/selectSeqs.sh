#! bin/bash
cp $0 input.fa;
grep ">" $1 | sed "s/>//g" > patterns.txt;
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' input.fa > tmp;
grep -Ff patterns.txt tmp > tmp2;
tr "\t" "\n" tmp2 > $0;
rm input.fa tmp tmp2;
