inputRedirect = {};
inputRedirect["01"] = "Universal";
inputRedirect["02"] = "/home/noor/Stokes_AntiStokes/RealPAP/BUSTED/RealPAPNumbered.txt";
inputRedirect["03"] = "/home/noor/Stokes_AntiStokes/RealPAP/BUSTED/RealPAPTree.trees";
inputRedirect["04"] = "All";
inputRedirect["05"] = "";

ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles/SelectionAnalyses/" + DIRECTORY_SEPARATOR + "BUSTED.bf", inputRedirect);
