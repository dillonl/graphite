"""
stolen from github.com/cc2qe/svtyper
"""
from math import log


def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for d in xrange(1, k + 1):
        r += log(n, 10)
        r -= log(d, 10)
        n -= 1
    return r

def gtl(ref, alt, is_dup, priors=(1e-4, 0.99, 0.01)):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [0.01, 0.3, 0.5]
    else:
        p_alt = [0.01, 0.5, 0.9]

    total = ref + alt

    lc_tot_alt = log_choose(total, alt)

    lp_homref = lc_tot_alt + alt * log(p_alt[0], 10) + ref * log(1 - p_alt[0], 10)#+ log(priors[0], 10)
    lp_het = lc_tot_alt + alt * log(p_alt[1], 10) + ref * log(1 - p_alt[1], 10)#+ log(priors[1], 10)
    lp_homalt = lc_tot_alt + alt * log(p_alt[2], 10) + ref * log(1 - p_alt[2], 10)#+ log(priors[2], 10)

    return (lp_homref, lp_het, lp_homalt)

def genotype(dp4):
    d = {}
    d['homref'], d['het'], d['homalt'] = gtl(sum(dp4[:2]), sum(dp4[2:]), False)
    gt, ll = 'homref', d['homref']
    for k in ('het', 'homalt'):
        if d[k] > ll:
            ll = d[k]
            gt = k

    s = dict(GT={"homref": "0/0", "het": "0/1", "homalt": "1/1"}[gt])
    s["GL"] = ",".join(("%.2f" % v for v in (d['homref'], d['het'], d['homalt'])))
    s["DP4"] = ",".join(map(str, dp4))
    s["DP"] = str(sum(dp4))
    return s

if __name__ == "__main__":

    import sys
    import toolshed as ts

    fh = ts.nopen(sys.argv[1])

    for header in ts.reader(fh, header=False):
        print "\t".join(header)
        if header[0] == "#CHROM": break

    header = header[:8]
    header.extend(["FORMAT", "sample"])
    fmt = "GT:GL:DP4:DP"
    for var in ts.reader(fh, header=header):
        var['FILTER'] = "."
        try:
            dp4 = [x for x in var['INFO'].split(";") if x.startswith("DP4=")][0][4:]
            dp4 = map(int, dp4.split(","))
        except IndexError:
            info = dict(x.split("=") for x in var['INFO'].split(";"))
            try:
                srf = int(info['SRF'])
                srr = int(info['SRR'])
                saf = sum(map(int, info['SAF'].split(",")))
                srf = sum(map(int, info['SRF'].split(",")))
                dp4 = [srf, srr, saf, srf]
            except:

                sinfo = dict(zip(var["FORMAT"].split(":"), var['sample'].split(":")))
                if "PL" in sinfo:
                    print "\t".join(var[h] for h in header)
                    continue

                ad = sinfo['AD']

                ad = ad[0], sum(ad[1:])
                assert len(ad) == 2, (ad, info, var)
                dp4 = [ad[0]/2, ad[0]/2, ad[1]/2, ad[1]/2]

                #print >>sys.stderr, var
                #print >>sys.stderr, info
                #print >>sys.stderr, info
                #raise
                #dp = int(info['DP'])
                #af = float(info['AF'])
                #dp4 = [0, int(dp * (1-af) + 0.5), 0, int(dp * af + 0.5)]

        samp = genotype(dp4)
        var["FORMAT"] = fmt
        var["sample"] = ":".join((samp["GT"], samp["GL"], samp["DP4"]))
        print "\t".join(var[h] for h in header)
