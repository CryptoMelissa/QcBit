def small_subset(n,k):
                ss = sage.combinat.subset.SubMultiset_sk(set(range(n)),k)
                return ss.random_element()
def hamming_weight(v):
                return sum((int(x) for x in v))

def circulant(v):
                l = len(v)
                return Matrix(GF(2),[[v[(i-j)%l] for i in xrange(l)] for j in xrange(l)])

def coeffs(p,sh=0):
    return ((p*x^sh).lift().padded_list() + [0]*n)[0:n]


def small(k):
    cs = set()
    for i in xrange(k):
        p = None
        while p is None or p in cs:
            p = randint(0,n-1)
        cs.add(p)
    return sum((x^i for i in cs))


def test(nn,w0,w1,block,tries=10000000):
                global gg
                global R
                global n
                global x
                global xr
                global RR
                n=nn
                R.<xr> = PolynomialRing(GF(2))
                RR.<x> = R.quotient(xr^n-1)
                mtries = 0
                print "Generating matrix..."
                while True:
                                try:
                                                h1 = small(w1)
                                                h0 = small(w0)
                                                g=h1/h0
                                                g = circulant(coeffs(g))
                                                h0 = Matrix(GF(2),coeffs(h0)).transpose()
                                                h1 = Matrix(GF(2),coeffs(h1)).transpose()
                                                break
                                except ZeroDivisionError:
                                                mtries = mtries + 1
                                                if mtries > 10:
                                                                raise Exception("Too many failures")


                

                assert g.transpose()*h0 == h1
                
                parity_equations=[]
                parity = []
                columns = []
                for b in xrange(0,n,block):
                                if any((b+j < n and h0[b+j] != 0 for j in xrange(block))):
                                                columns.extend(range(b,min(n,b+block)))
                                                v = [0]*n
                                                par=0
                                                for j in xrange(block) :
                                                    if b+j <n :
                                                        v[b+j] = 1
                                                    par=par+h0[b+j][0]
                                                parity.append(par)
                                                parity_equations.append(vector(v))
                                                
                nb_parity_equations = len(parity_equations)
                equations = g.transpose().rows()
                equations.extend(parity_equations)
                G = Matrix(GF(2),equations).transpose()
                
                
                
                print "   Matrix is generated."
                print "Expected number of tries:", N(exp(w1*(len(columns)-nb_parity_equations)/n))
                
                
                G = G[columns].transpose()
                g=g[columns].transpose()
                columns = set(columns)
                right_member = Matrix(GF(2),[0]*(len(columns)-nb_parity_equations) + parity).transpose()
                for i in xrange(tries):
                                #print "Try #%d" % i
                                subset = small_subset(n,len(columns)-nb_parity_equations)
                                subset.extend(range(n,(n+nb_parity_equations)))
                                gg = G[subset]
                                try:
                                    for v in gg.transpose().solve_left(right_member.transpose()):
                                                h0r = []
                                                if hamming_weight(v) == w0 and hamming_weight(g*v) == w1:
                                                                print "Found result at test number %d" %i
                                                                off = 0
                                                                for j in xrange(n):
                                                                                if j in columns:
                                                                                                h0r.append(v[off])
                                                                                                off += 1
                                                                                else: h0r.append(0)
                                                                h0r = Matrix(GF(2),[h0r]).transpose()
                                                                h1r = Matrix(g*v).transpose()
                                                                
                                                                if h0r == h0: print "  h0 is correct"
                                                                if h1r == h1: print "  h1 is correct"														 if (h0r == h0) and (h1r == h1): return
                                except ValueError:
                                    continue