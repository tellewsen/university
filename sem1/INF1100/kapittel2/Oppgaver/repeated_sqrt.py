from math import sqrt
for n in range(1,60):
    r = 2.0
    for i in range(n):
        r = sqrt(r)
    for i in range(n):
        r=r**2
    print '%d times sqrt and **2: %.16f' % (n,r)

"""
Terminal> Oppgaver $ python repeated_sqrt.py 
1 times sqrt and **2: 2.0000000000000004
2 times sqrt and **2: 1.9999999999999996
3 times sqrt and **2: 1.9999999999999996
4 times sqrt and **2: 1.9999999999999964
5 times sqrt and **2: 1.9999999999999964
6 times sqrt and **2: 1.9999999999999964
7 times sqrt and **2: 1.9999999999999714
8 times sqrt and **2: 2.0000000000000235
9 times sqrt and **2: 2.0000000000000235
10 times sqrt and **2: 2.0000000000000235
11 times sqrt and **2: 2.0000000000000235
12 times sqrt and **2: 1.9999999999991336
13 times sqrt and **2: 1.9999999999973292
14 times sqrt and **2: 1.9999999999973292
15 times sqrt and **2: 1.9999999999973292
16 times sqrt and **2: 2.0000000000117746
17 times sqrt and **2: 2.0000000000408580
18 times sqrt and **2: 2.0000000000408580
19 times sqrt and **2: 2.0000000001573586
20 times sqrt and **2: 2.0000000001573586
21 times sqrt and **2: 2.0000000001573586
22 times sqrt and **2: 2.0000000010885857
23 times sqrt and **2: 2.0000000029511749
24 times sqrt and **2: 2.0000000066771721
25 times sqrt and **2: 2.0000000066771721
26 times sqrt and **2: 1.9999999917775542
27 times sqrt and **2: 1.9999999917775542
28 times sqrt and **2: 1.9999999917775542
29 times sqrt and **2: 1.9999999917775542
30 times sqrt and **2: 1.9999999917775542
31 times sqrt and **2: 1.9999999917775542
32 times sqrt and **2: 1.9999990380770896
33 times sqrt and **2: 1.9999971307544144
34 times sqrt and **2: 1.9999971307544144
35 times sqrt and **2: 1.9999971307544144
36 times sqrt and **2: 1.9999971307544144
37 times sqrt and **2: 1.9999971307544144
38 times sqrt and **2: 1.9999360966436217
39 times sqrt and **2: 1.9999360966436217
40 times sqrt and **2: 1.9999360966436217
41 times sqrt and **2: 1.9994478907329654
42 times sqrt and **2: 1.9984718365144798
43 times sqrt and **2: 1.9965211562778555
44 times sqrt and **2: 1.9965211562778555
45 times sqrt and **2: 1.9887374575497223
46 times sqrt and **2: 1.9887374575497223
47 times sqrt and **2: 1.9887374575497223
48 times sqrt and **2: 1.9887374575497223
49 times sqrt and **2: 1.8682459487159784
50 times sqrt and **2: 1.6487212645509468
51 times sqrt and **2: 1.6487212645509468
52 times sqrt and **2: 1.0000000000000000
53 times sqrt and **2: 1.0000000000000000
54 times sqrt and **2: 1.0000000000000000
55 times sqrt and **2: 1.0000000000000000
56 times sqrt and **2: 1.0000000000000000
57 times sqrt and **2: 1.0000000000000000
58 times sqrt and **2: 1.0000000000000000
59 times sqrt and **2: 1.0000000000000000
"""


"""
Program takes sqrt of 2 and then the answer of that **2
Because of the way a computer stores numbers we don't get the exact answer. 
This error gets worse by each time we do this and after 52 runs, the answer is completely wrong. It gets bad long before that.
"""
