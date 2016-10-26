           REAL FUNCTION FUNC(X)
           COMMON/PAWPAR/PAR(3)
           FUNC=PAR(1)*((PAR(2)**(x/PAR(3)))/(Gamma(x/PAR(3)+1)))*
     &	   exp(-1.*PAR(2))
           END
