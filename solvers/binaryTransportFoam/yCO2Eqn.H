{
   // surfaceScalarField PHf = fvc::interpolate(PH,"PHf");
    fvScalarMatrix yCO2Eqn
        (
            eps*(fvm::ddt(yCO2) - fvm::ddt(Sb,yCO2)) + fvm::div(phia,yCO2)
            ==
            // event source terms
            - sourceTerm
        );
    yCO2Eqn.solve();
    yH2O= 1-yCO2;
    
    //NG = eps*rhoa*(1-Sb)/(Mg); NCO2g = yCO2*NG; NG = NCO2g; //NG_new = NCO2g = eps*rhoa*(1-Sb_new)/Mg = yCO2*eps*rhoa*(1-Sb_old)/(Mg); => Sg = yCO2*Sg => Sw = 1-yCO2*Sg = 1- yCO2*(1-Sw_old) = 1-yCO2+yCO2*Sw_old;
    //forAll(yCO2, cellI)
    //{
	//    if (yCO2[cellI]> 1)//Sb.time()[0][cellI]) //Sb.internalField()[cellI]) //Sb[0][cellI]) //Sb.time().timeIndex(0)) //Sb[0]) //Sb.time().timeIndex()) //Sb.time().valueList()[0]) //Sb.time().index(0))//(xCO2[cellI]==0)
	//	    {
	//		    Sb[cellI] = 1-yCO2[cellI]*(1-Sb[cellI]); //Sg = yCO2*Sg; Sg = 1-Sw = yCO2*Sg = yCO2*(1-Sw_old);
	//		    yCO2[cellI] =1;
	//		    yH2O[cellI]= 0;
	//	    }
	//	    else if (yCO2[cellI]< 0)////NG = eps*rhoa*(1-Sb)/(Mg); NH2Og = yH2O*NG; NG = NH2Og; //NG_new = NH2Og = eps*rhoa*(1-Sb_new)/Mg = yH2O*eps*rhoa*(1-Sb_old)/(Mg); => Sg = yH2O*Sg => Sw = 1-yH2O*Sg = 1- yH2O*(1-Sw_old) = 1-yH2O+yH2O*Sw_old;
	//	    {
	//		    Sb[cellI] = 1-yH2O[cellI]*(1-Sb[cellI]); //Sg = yCO2*Sg; Sg = 1-Sw = yCO2*Sg = yCO2*(1-Sw_old);
	//		    yH2O[cellI] =1;
	//		    yCO2[cellI]= 0;
	//	    }
    //}
     

    Info << "yH2O: " << " Min(yH2O) = " << gMin(yH2O) << " Max(yH2O) = " << gMax(yH2O) << endl;

}
