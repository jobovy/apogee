{

count = count + 1

if (count == 1) print "KURUCZ"
if (count == 1) printf("%s", $0)
if (count == 2) printf("%s\n", $0)

if ($3 == "GRAVITY") logg = $4

if ($1 == "ABUNDANCE" && $2 == "SCALE") zero = log($3)/log(10.0)

if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $9 == "6" ) c = 12.0 + $10 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $11 == "7" ) n = 12.0 + $12 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $13 == "8" ) o = 12.0 + $14 + zero + 0.04

if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $7 == "11" ) na = 12.0 + $8 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $9 == "12" ) mg = 12.0 + $10 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $11 == "13" ) al  = 12.0 + $12 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $13 == "14" ) si = 12.0 + $14 + zero + 0.04

if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $5 == "16" ) s  = 12.0 + $6 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $11 == "19" ) k  = 12.0 + $12 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $13 == "20" ) ca = 12.0 + $14 + zero + 0.04

if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $5 == "22" ) ti  = 12.0 + $6 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $7 == "23" ) v = 12.0 + $8 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $9 == "24" ) cr = 12.0 + $10 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $11 == "25" ) mn  = 12.0 + $12 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $13 == "26" ) fe = 12.0 + $14 + zero + 0.04

if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $3 == "27" ) co  = 12.0 + $4 + zero + 0.04
if ($1 == "ABUNDANCE" && $2 == "CHANGE" && $5 == "28" ) ni  = 12.0 + $6 + zero + 0.04

if (start > 0 && start <= lines)  {
    print $0
    start = start +1
    micro = $7 * 1.0
    }

if ($1 == "READ") lines = $3
if ($1 == "READ") start = 1
if ($1 == "READ") printf("             %s\n",$3)

#micro = micro   # substitute your favorite formula here
micro = (2.24 - 0.3 * logg) * 1E5   # substitute your favorite formula here
if (length(vmicro)>0) micro=vmicro * 1E5

}
END{printf("%10.3f\n",micro / 1E5)
printf("NATOMS     17 %6.3f\n",zero)
printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",6.0,c,7.0,n,8.0,o)
printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",11.0,na,12.0,mg,13.0,al,14.0,si)
printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",16.0,s,19.0,k,20.0,ca)
printf("%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f%10.1f%10.3f\n",22.0,ti,23.0,v,24.0,cr,25.0,mn,26.0,fe)
printf("%10.1f%10.3f%10.1f%10.3f\n",27.0,co,28.0,ni)
printf("NMOL          20\n")
printf("       1.1     107.0     108.0     607.0     608.0     708.0       6.1\n")
printf("       7.1       8.1      12.1     112.0     101.0     106.0     101.0\n")
printf("      22.1     822.0      14.1     114.0      26.1     126.0\n")
} 
