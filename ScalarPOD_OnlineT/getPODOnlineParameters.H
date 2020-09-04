struct PODOnlineParameters
{
    List<fileName> folders_list;
    word fieldName;
    label  BasisNumber;
    bool isMultiDir;
};

void readfolderList
(
    const IOdictionary& PODsolverDict,
    List<fileName>& folders_list
)
{
    List<fileName> tmpfolderList (PODsolverDict.subDict("Online_parameters").lookup("foldersList"));
    folders_list = tmpfolderList;
};

PODOnlineParameters getPODOnlineParameters (Foam::argList args)
{

    Foam::Time runTime(Foam::Time::controlDictName, args);
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );


    IOdictionary PODsolverDict
    (
        IOobject
        (
            "PODsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read Initial and last time from the POD dictionary
    const entry* existField = PODsolverDict.subDict("Online_parameters").lookupEntryPtr("field", false, true);
    const entry* existfoldersLists = PODsolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existNumberBasis = PODsolverDict.subDict("Online_parameters").lookupEntryPtr("BasisNumber", false, true);

    // Initiate variable from PODsolverDict
    if (!(existNumberBasis) )
    {
        Info << "Error: you must define number of Basis: 'BasisNumber'" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' " << endl;
        abort();
    }


    List<fileName> folders_list;

    bool isMultiDir;

    if (!existfoldersLists)
    {
        folders_list.setSize(1);
        folders_list[0]=PODsolverDict.caseName();
        Info << "Reading list of folders\n" << folders_list << endl;
        isMultiDir=false;
    }
    else
    {
        readfolderList(PODsolverDict,folders_list);
        Info << "Reading list of folders\n" << folders_list << endl;
        isMultiDir=true;
    }

    word fieldName
    (
        PODsolverDict.subDict("Online_parameters").lookup("field")
    );

    label BasisNumber =
        readInt
        (
            PODsolverDict.subDict("Online_parameters").lookup("BasisNumber")
        );
    Info << "Number of Basis: " << BasisNumber << endl;



    PODOnlineParameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.BasisNumber=BasisNumber;
    results.isMultiDir=isMultiDir;

    return results;
}

