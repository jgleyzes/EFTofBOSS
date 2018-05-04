#include "RedshiftBiasEFT.h"

static YesNo GetAnswer (const string & Answer) {
	if (Answer == "Y" || Answer == "y" || Answer == "yes" || Answer == "Yes" || Answer == "YES") return 1 ;
	if (Answer == "N" || Answer == "n" || Answer == "no" || Answer == "No" || Answer == "NO") return 0 ;
	else return 1 ; 
}

void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToFolder, string & PathToFolderRD,string & PathToFolderCosmoRef, string & PathToLinearPowerSpectrum, string & PathToTriangles, 
	YesNo & ComputePowerSpectrum, YesNo & ComputeBispectrum, YesNo & UseRef, YesNo & ImportM, YesNo & ExportM, 
	PrecisionIntegr & Eps, double & EpsRel_IntegrBispectrumAP, double & aperp, double & apar) {

///// Reading ini file /////
	ifstream inFile (ConfigFile) ;    
	char line[1024] ;
	char delim[] = " =" ;
	char * token ;

	bool physical_b = false, physical_cdm = false ;
	
	if (!inFile.is_open()) {
		cout << "Could not open the input file: " << ConfigFile << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	while (inFile) {
		
		inFile.getline(line, 1024) ;
		
		if (strlen(line) == 0 or line[0]=='#' or line[0] == '[') continue ;
		
		token = strtok(line, delim) ;
		
		if (strcmp(token, "PathToFolder") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToFolder = token ;
		else if (strcmp(token, "PathToFolderRD") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToFolderRD = token ;
		else if (strcmp(token, "PathToFolderCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToFolderCosmoRef = token ;
		else if (strcmp(token, "PathToLinearPowerSpectrum") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToLinearPowerSpectrum = token ;
		else if (strcmp(token, "PathToTriangles") == 0) while ((token = strtok(NULL, delim)) != NULL) PathToTriangles = token ;

		else if (strcmp(token, "knl") == 0) while ((token = strtok(NULL, delim)) != NULL) knl = atof(token) ;
		else if (strcmp(token, "km") == 0) while ((token = strtok(NULL, delim)) != NULL) km = atof(token) ;
		else if (strcmp(token, "nbar") == 0) while ((token = strtok(NULL, delim)) != NULL) nbar = atof(token) ;
		else if (strcmp(token, "ComputePowerSpectrum") == 0) while ((token = strtok(NULL, delim)) != NULL) ComputePowerSpectrum = GetAnswer(token) ;
		else if (strcmp(token, "ComputeBispectrum") == 0) while ((token = strtok(NULL, delim)) != NULL) ComputeBispectrum = GetAnswer(token) ;
		else if (strcmp(token, "UseCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) UseRef = GetAnswer(token) ;
		else if (strcmp(token, "ImportResummationMatrix") == 0) while ((token = strtok(NULL, delim)) != NULL) ImportM = GetAnswer(token) ;
		else if (strcmp(token, "ExportResummationMatrix") == 0) while ((token = strtok(NULL, delim)) != NULL) ExportM = GetAnswer(token) ;
		
		else if (strcmp(token, "z_pk") == 0) while ((token = strtok(NULL, delim)) != NULL) z0 = atof(token) ;
		else if (strcmp(token, "ln10^{10}A_s") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[0] = atof(token) ;
		else if (strcmp(token, "n_s") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[1] = atof(token) ;
		else if (strcmp(token, "h") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[2] = atof(token) ;
		else if (strcmp(token, "Omega_b") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[3] = atof(token) ;
		else if (strcmp(token, "Omega_cdm") == 0) while ((token = strtok(NULL, delim)) != NULL) cosmo[4] = atof(token) ;
		else if (strcmp(token, "omega_b") == 0) while ((token = strtok(NULL, delim)) != NULL) { physical_b = true ; cosmo[3] = atof(token) ; }
		else if (strcmp(token, "omega_cdm") == 0) while ((token = strtok(NULL, delim)) != NULL) { physical_cdm = true ; cosmo[4] = atof(token) ; }

		else if (strcmp(token, "EpsAbs_NoCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) Eps[0][0] = atof(token) ;
		else if (strcmp(token, "EpsRel_NoCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) Eps[0][1] = atof(token) ;
		else if (strcmp(token, "EpsAbs_YesCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) Eps[1][0] = atof(token) ;
		else if (strcmp(token, "EpsRel_YesCosmoRef") == 0) while ((token = strtok(NULL, delim)) != NULL) Eps[1][1] = atof(token) ;

		else if (strcmp(token, "EpsRel_IntegrBispectrumAP") == 0) while ((token = strtok(NULL, delim)) != NULL) EpsRel_IntegrBispectrumAP = atof(token) ;
		else if (strcmp(token, "aperp") == 0) while ((token = strtok(NULL, delim)) != NULL) aperp = atof(token) ;
		else if (strcmp(token, "apar") == 0) while ((token = strtok(NULL, delim)) != NULL) apar = atof(token) ;



		}
	
	inFile.close() ;

	// Taking care of both fractional and physical densities
	if (physical_b == true) cosmo[3]/pow(cosmo[2],2) ;
	if (physical_cdm == true) cosmo[4]/pow(cosmo[2],2) ;
}

