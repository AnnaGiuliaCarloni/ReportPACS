/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"    	
#include "simpleControl.H"

#include "fvOptionList.H"  	

#include "IOdictionary.H"
#include "autoPtr.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "fvOptions.H"  
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            T.correctBoundaryConditions();{
            fvScalarMatrix TEqn
            (
                Foam::pow(L.value(),2) * fvm::ddt(T) - fvm::laplacian(DT, T) == frate*heatSource * Foam::pow(L.value(),2)
            );

            //fvOptions.constrain(TEqn);
            TEqn.solve();
            //fvOptions.correct(T);
            }
           
        }

	Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
	    // Dentro agli include ci sono risolte le equazioni che mi servono?
	
        //{
         
	#include "cEqn.H"
            
        //}


    runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
	
	
    //Info<< "Average Concentration = " << C.weightedAverage(mesh.V()) << endl;

    //dimensionedScalar averagedValue=C.weightedAverage(mesh.V());
    //Info<< "Average Concentration = " <<averagedValue.value() << endl;
 
    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
