struct PODparameters
{
    List<fileName> folders_list; //filename dovrebbe essere il tipo non ancora definito
    label  maxBasis;
    scalar accuracy;
    bool isMultiDir;
};

void readfolderList  
(
    const IOdictionary& PODsolverDict,
    List<fileName>& folders_list
)
{
    List<fileName> tmpfolderList (PODsolverDict.subDict("Offline_parameters").lookup("foldersList"));
    folders_list = tmpfolderList;
};

PODparameters getPODparameters (Foam::argList args)
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
    const entry* existfoldersLists = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existmaxBasis = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("maxBasis", false, true);
    const entry* existAccuracy = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("accuracy", false, true);

    // Initiate variable from PODsolverDict
    if (!(existmaxBasis) | !(existAccuracy))
    {
        Info << "Error: you must define accuracy and max number of Basis" << endl;
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
        PODsolverDict.subDict("Offline_parameters").lookup("field")
    );

    label maxBasis =
        readInt
        (
            PODsolverDict.subDict("Offline_parameters").lookup("maxBasis")
        );
    Info << "Seeking number of Basis: " << maxBasis << endl;

    scalar accuracy =
        readScalar
        (
            PODsolverDict.subDict("Offline_parameters").lookup("accuracy")
        );

    Info << "Seeking accuracy: " << accuracy << endl;



    PODparameters results;
    results.folders_list=folders_list;
    results.maxBasis=maxBasis;
    results.accuracy=accuracy;
    results.isMultiDir=isMultiDir;

    return results;
}

