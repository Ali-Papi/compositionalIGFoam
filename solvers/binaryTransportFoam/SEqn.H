{

    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();
    forAll(mesh.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Ua.boundaryField()[patchi]))
        {
            phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Ub.boundaryField()[patchi]))
        {
            phib.boundaryFieldRef()[patchi] = Ub.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
  
    fvScalarMatrix SbEqn
        (
            -eps*fvm::ddt(Sb) + fvc::div(phia) 
            ==
            // event source terms
            - sourceTerm
        );

    SbEqn.solve();
    
    //forAll(Sb, cellI)
    //{
	//    if (Sb[cellI] <Sbmin1)//Sb.time()[0][cellI]) //Sb.internalField()[cellI]) //Sb[0][cellI]) //Sb.time().timeIndex(0)) //Sb[0]) //Sb.time().timeIndex()) //Sb.time().valueList()[0]) //Sb.time().index(0))//(xCO2[cellI]==0)
	//	    {
	//		    Info << "Sbmin: " << Sbmin1 << endl;
	//		    Sb[cellI] = 0.1; //2*Sbmin1; //Sg = yCO2*Sg; Sg = 1-Sw = yCO2*Sg = yCO2*(1-Sw_old);
	//	    }
   // }
    
    Info << "Saturation b: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
