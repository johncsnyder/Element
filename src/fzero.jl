

doc"""
`fzero(f,x1,x2,args...;tol=1e-12,maxsteps=100)`

finds the root $x_0$ of $f(x)$ in the interval $[x_1,x_2]$, such that
$f(x_0) = 0$ (up to specified tolerance `tol`).

`f` - a function f(x,args...) that returns Float64, where extra parameters 
may be passed through to `f`.

throws an error if $f(x_1)$ and $f(x_2)$ have the same sign, or 
if more than `maxsteps` iterations are taken.

uses the Dekker-Brent method (see `Numerical Recipes in C`, p. 361).
"""
@generated function fzero(f, x1::Float64, x2::Float64, args...; tol::Float64=1e-12, maxsteps::Int=100)
    quote
        local i,d,e,min1,min2,fc,p,q,r,s,tol1,xm
        a = x1
        b = x2
        c = x2
        fa = f(a,args...)::Float64
        fb = f(b,args...)::Float64
        δ = 3e-12
        ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) && error("f(a) and f(b) must have opposite signs")
        fc = fb
        for i in 1:maxsteps
            if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)
                c = a
                fc = fa
                e = d = b-a
            end
            if abs(fc) < abs(fb)
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end
            tol1 = 2*δ*abs(b) + 0.5*tol
            xm = 0.5*(c-b)
            if abs(xm) <= tol1 || fb == 0.0
                return b
            end
            if abs(e) >= tol1 && abs(fa) > abs(fb)
                s = fb/fa
                if a == c
                    p = 2.0*xm*s
                    q = 1.0-s
                else
                    q = fa/fc
                    r = fb/fc
                    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                    q = (q-1.0)*(r-1.0)*(s-1.0)
                end
                if p > 0.0
                    q = -q
                end
                p = abs(p)
                min1 = 3.0*xm*q - abs(tol1*q)
                min2 = abs(e*q)
                if 2.0*p < (min1 < min2 ? min1 : min2)
                    e = d
                    d = p/q
                else
                    d = xm
                    e = d
                end
            else
                d = xm
                e = d
            end
            a = b
            fa = fb
            if abs(d) > tol1
                b += d
            else
                b += abs(tol1)*sign(xm)
            end
            fb = f(b,args...)::Float64
        end
        error("maximum number of iterations exceeded")
    end
end








