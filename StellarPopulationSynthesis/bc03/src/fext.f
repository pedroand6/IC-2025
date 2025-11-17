        REAL FUNCTION FEXT(X,TV)

	! Attenuation function used in Charlot and Fall (2000) model

        real x,tv,tau

        fext=1.

        if (tv.eq.0.) return

        tau=tv*( (5500./x)**0.7 )
        fext=exp(-tau)

        return
        end
