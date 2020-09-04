/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    ScalarPOD_Online_Thesis

Description

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "MOR.H"
#include "simpleMatrix.H"
#include "IOmanip.H"
#include "IOstream.H"


// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "getPODOnlineParameters.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "addDictOption.H"
    timeSelector::addOptions(false, true);

#include "setRootCase.H"
Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"
    //get POD_Online parameters

    PODOnlineParameters POD_parameters = getPODOnlineParameters(args);

   
    //List<fileName> foldersList = POD_parameters.folders_list ;

    label BasisNumber= POD_parameters.BasisNumber;


    /***********************************************************************************/

    PtrList<volScalarField>  scalarPODOrthoNormalBasis (BasisNumber);

    bool PODBasisLoaded = false;

    while (PODBasisLoaded == false)
    {
        //Foam::Time runTime(args.rootPath(), args.caseName());
//#include "CreateMesh.H"
        //caricamento basi
        
        runTime.setTime(0, 0);
        for(label basisI=0; basisI< BasisNumber; basisI ++)
        {
            //runTime.setTime(mfI, mfI);
            scalarPODOrthoNormalBasis.set
            (
                basisI,
                new volScalarField
                (
                    IOobject
                    (
                        "CPOD"+name(basisI),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        PODBasisLoaded = true;
    }

//caricamento dei coefficienti

scalarRectangularMatrix COEFF(IFstream(args.rootPath()/args.caseName()/"coefficientiPOD/coeffs_out1.txt")()); 

scalarRectangularMatrix TIMESTEPS(IFstream(args.rootPath()/args.caseName()/"coefficientiPOD/t_out1.txt")());

//Info << TIMESTEPS << endl;



    /***********************************************************************************/

//Ricostruzione delle soluzioni(ultima modifica)

for (label timeIndex = 0; timeIndex<TIMESTEPS.n(); ++timeIndex)
{

    volScalarField reconstructed
    (
        IOobject
        (
            "CPODreconstruct"+name(TIMESTEPS[0][timeIndex]), //name trasforma numero in stringa
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensioned<scalar>
        (
            "zero",
            scalarPODOrthoNormalBasis[0].dimensions(), //le dimensioni delcampo ricostruito sono quelle della T
            pTraits<scalar>::zero
        )
    );

    for (label baseI = 0; baseI < BasisNumber; baseI++)
    {
        reconstructed +=
            COEFF[timeIndex][baseI]*scalarPODOrthoNormalBasis[baseI];
    }

reconstructed.write();

}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End" << endl;

}


// ************************************************************************* //
