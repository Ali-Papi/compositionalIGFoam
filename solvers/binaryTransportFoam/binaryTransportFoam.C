/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    impesFoam

Description
    Transient solver for incompressible two-phase flow (Darcy's law) in porous media
    using the IMPES method (IMplicit Pressure Explicit Saturation).
    Permeability is isotropic (K == volScalarField)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "harmonic.H"
#include "incompressiblePhase.H"
#include "twophasePorousMediumModel.H"
#include "capillarityModel.H"
#include "relativePermeabilityModel.H"
#include "sourceEventFile.H"
#include "outputEventFile.H"
#include "patchEventFile.H"
#include "eventInfiltration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "../headerPMF.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createSbFields.H"
    #include "readTimeControls.H"
    #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (sourceEventIsPresent) sourceEvent.updateIndex(runTime.userTimeValue());
        forAll(patchEventList,patchEventi) patchEventList[patchEventi]->updateIndex(runTime.userTimeValue());
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "computeSourceTerm.H"
				
				//SbLastItr1 = Sb;
        //Sb_error = 1;
				//while (Sb_error>1e-5)
				//{
			//- Solve saturation equation (explicit)
					#include "SEqn.H"
					#include "updateSbProperties.H"

				#include "yCO2Eqn.H"
				#include "xCO2Eqn.H"
          //- Solve pressure equation (implicit)
          #include "pEqn.H"
				
				  //Sb_error= max(Sb.field()-SbLastItr.field());
					//SbLastItr = Sb;
					//Info << "Sb_error is = " << Sb_error << endl;
				//}
				
				//?? #include "updateSbProperties.H"
				

        #include "eventWrite.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
