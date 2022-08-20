for n in range(0,80):
    a = 1.0 +2.0**-n
    if a==1:
        print 'n=%g   a=%.16f' %(n,a)
        print '1 + 2**-%g gir feil resultat' %n
        break
