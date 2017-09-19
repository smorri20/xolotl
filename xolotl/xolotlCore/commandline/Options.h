#ifndef OPTIONS_H
#define OPTIONS_H

#include "IOptions.h"
#include "optionhandlers/IOptionHandler.h"

namespace xolotlCore {

/**
 * Options realizes the IOptions interface.
 * All private members will be accessed through getters and setters.
 */
class Options: public IOptions {

protected:
	/**
	 * Map of options we support, keyed by option switch string.
	 */
	typedef std::map<std::string, IOptionHandler*> OptionsMap;
	OptionsMap optionsMap;

	/**
	 * The flag that says if Xolotl should run.
	 */
	bool shouldRunFlag;

	/**
	 * The value of the exit code. Should be 0 if everything went well.
	 */
	int exitCode;

	/**
	 * The name of the file where the network is stored.
	 */
	std::string networkFilename;

	/**
	 * The number of options that will be given to PETSc.
	 */
	int petscArgc;

	/**
	 * The pointer to the options that will be given to PETSc.
	 */
	char **petscArgv;

	/**
	 * Use the constant temperature set of handlers?
	 */
	bool constTempFlag;

	/**
	 * Value for the constant temperature.
	 */
	double constTemperature;

	/**
	 * Value for the temperature gradient
	 */
	double temperatureGradient;

	/**
	 * Use the temperature profile set of handlers?
	 */
	bool tempProfileFlag;

	/**
	 * Name of the input temperature profile file.
	 */
	std::string tempProfileFilename;

	/**
	 * Use the flux amplitude option?
	 */
	bool fluxFlag;

	/**
	 * Value for the  flux.
	 */
	double fluxAmplitude;

	/**
	 * Use a time profile for the flux?
	 */
	bool fluxProfileFlag;

	/**
	 * Name of the input time profile file for the flux.
	 */
	std::string fluxProfileFilename;

	/**
	 * Which type of performance infrastructure should we use?
	 */
	xolotlPerf::IHandlerRegistry::RegistryType perfRegistryType;

	/**
	 * Use the "standard" set of handlers for the visualization infrastructure?
	 */
	bool vizStandardHandlersFlag;

	/**
	 * Name of the material.
	 */
	std::string materialName;

	/**
	 * Value of the initial vacancy concentration.
	 */
	double initialVConcentration;

	/**
	 * Value of the portion of the void on the grid at the start of the simulation.
	 */
	double voidPortion;

	/**
	 * Number of dimensions for the simulation.
	 */
	int dimensionNumber;

	/**
	 * Use a regular grid on the x direction?
	 */
	bool useRegularGridFlag;

	/**
	 * The map of physical processes to use in the simulation.
	 */
	std::map<std::string, bool> processMap;

	/**
	 * String of the list of wanted GB.
	 */
	std::string gbList;

	/**
	 * Minimum size for the grouping.
	 */
	int groupingMin;

	/**
	 * Width for the grouping in the first direction.
	 */
	int groupingWidthA;

	/**
	 * Width for the grouping in the second direction.
	 */
	int groupingWidthB;

	/**
	 * Value of the sputtering yield.
	 */
	double sputteringYield;

	/**
	 * Use a HDF5 file?
	 */
	bool useHDF5Flag;

	/**
	 * Use the phase cut for the network?
	 */
	bool usePhaseCutFlag;

	/**
	 * Maximum number of He or Xe
	 */
	int maxImpurity;

	/**
	 * Maximum number of V
	 */
	int maxV;

	/**
	 * Maximum number of I
	 */
	int maxI;

	/**
	 * Number of grid point in the depth direction
	 */
	int nX;

	/**
	 * Step size in the depth direction
	 */
	double xStepSize;

	/**
	 * Number of grid point in the Y direction
	 */
	int nY;

	/**
	 * Step size in the Y direction
	 */
	double yStepSize;

	/**
	 * Number of grid point in the Z direction
	 */
	int nZ;

	/**
	 * Step size in the Z direction
	 */
	double zStepSize;

    /**
     * Whether to write reaction network to debug output file.
     */
    bool shouldWriteDebugNetwork;

    /**
     * Name of output file for debug reaction network.
     * Ignored unless shouldWriteDebugNetwork is true.
     */
    std::string debugNetworkFilename;

public:

	/**
	 * The constructor.
	 */
	Options();

	/**
	 * The destructor.
	 */
	~Options();

	/**
	 * Read the parameters from the given file to set the different
	 * xolotl options.
	 * \see IOptions.h
	 */
	void readParams(char* argv[]) override;

	/**
	 * Show our help message.
	 * \see IOptions.h
	 */
	void showHelp(std::ostream& os) const override;

	/**
	 * Should the program run after parsing the parameter file?
	 * \see IOptions.h
	 */
	bool shouldRun() const override {
		return shouldRunFlag;
	}

	/**
	 * Set the shouldRunFlag.
	 * \see IOptions.h
	 */
	void setShouldRunFlag(bool flag) override {
		shouldRunFlag = flag;
	}

	/**
	 * If program shouldn't run, what should its exit code be?
	 * \see IOptions.h
	 */
	int getExitCode() const override {
		return exitCode;
	}

	/**
	 * Set the value for the exit code.
	 * \see IOptions.h
	 */
	void setExitCode(int code) override {
		exitCode = code;
	}

	/**
	 * Get the name of the network file.
	 * \see IOptions.h
	 */
	std::string getNetworkFilename() const override {
		return networkFilename;
	}

	/**
	 * Set the name of the network file.
	 * \see IOptions.h
	 */
	void setNetworkFilename(const std::string& name) override {
		networkFilename = name;
	}

	/**
	 * Get the Argc for PETSc.
	 * \see IOptions.h
	 */
	int getPetscArgc() const override {
		return petscArgc;
	}

	/**
	 * Set the Argc for PETSc.
	 * \see IOptions.h
	 */
	void setPetscArgc(int argc) override {
		petscArgc = argc;
	}

	/**
	 * Get the Argv for PETSc.
	 * \see IOptions.h
	 */
	char** getPetscArgv() const override {
		return petscArgv;
	}

	/**
	 * Set the Argv for PETSc.
	 * \see IOptions.h
	 */
	void setPetscArgv(char** argv) override {
		petscArgv = argv;
	}

	/**
	 * Should we use const temperature handlers?
	 * \see IOptions.h
	 */
	bool useConstTemperatureHandlers() const override {
		return constTempFlag;
	}

	/**
	 * Set the constTempFlag.
	 * \see IOptions.h
	 */
	void setConstTempFlag(bool flag) override {
		constTempFlag = flag;
	}

	/**
	 * Obtain the value of the constant temperature to be used.
	 * \see IOptions.h
	 */
	double getConstTemperature() const override {
		return constTemperature;
	}

	/**
	 * Set the constant temperature.
	 * \see IOptions.h
	 */
	void setConstTemperature(double temp) override {
		constTemperature = temp;
	}

	/**
	 * Obtain the value of the temperature gradient to be used.
	 * \see IOptions.h
	 */
	double getTemperatureGradient() const override {
		return temperatureGradient;
	}

	/**
	 * Set the temperature gradient.
	 * \see IOptions.h
	 */
	void setTemperatureGradient(double grad) override {
		temperatureGradient = grad;
	}

	/**
	 * Should we use temperature profile handlers?
	 * \see IOptions.h
	 */
	bool useTemperatureProfileHandlers() const override {
		return tempProfileFlag;
	}

	/**
	 * Set the tempProfileFlag.
	 * \see IOptions.h
	 */
	void setTempProfileFlag(bool flag) override {
		tempProfileFlag = flag;
	}

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 * \see IOptions.h
	 */
	std::string getTempProfileFilename() const override {
		return tempProfileFilename;
	}

	/**
	 * Set the name of the profile file to use.
	 * \see IOptions.h
	 */
	void setTempProfileFilename(const std::string& name) override {
		tempProfileFilename = name;
	}

	/**
	 * Should we use the flux option?
	 * \see IOptions.h
	 */
	bool useFluxAmplitude() const override {
		return fluxFlag;
	}
	;

	/**
	 * Set the fluxFlag.
	 * \see IOptions.h
	 */
	void setFluxFlag(bool flag) override {
		fluxFlag = flag;
	}

	/**
	 * Obtain the value of the flux intensity to be used.
	 * \see IOptions.h
	 */
	double getFluxAmplitude() const override {
		return fluxAmplitude;
	}

	/**
	 * Set the value for the flux intensity to use.
	 * \see IOptions.h
	 */
	void setFluxAmplitude(double flux) override {
		fluxAmplitude = flux;
	}

	/**
	 * Should we use a time profile for the flux?
	 * \see IOptions.h
	 */
	bool useFluxTimeProfile() const override {
		return fluxProfileFlag;
	}

	/**
	 * Set the fluxProfileFlag.
	 * \see IOptions.h
	 */
	void setFluxProfileFlag(bool flag) override {
		fluxProfileFlag = flag;
	}

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * flux.
	 * \see IOptions.h
	 */
	std::string getFluxProfileName() const override {
		return fluxProfileFilename;
	}

	/**
	 * Set the name of the time profile file to use.
	 * \see IOptions.h
	 */
	void setFluxProfileName(const std::string& name) override {
		fluxProfileFilename = name;
	}

	/**
	 * Which type of performance handlers should we use?
	 * \see IOptions.h
	 */
	xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(void) const override {
		return perfRegistryType;
	}

	/**
	 * Set the type of performance handlers to use.
	 * \see IOptions.h
	 */
	void setPerfHandlerType(xolotlPerf::IHandlerRegistry::RegistryType rtype) override {
		perfRegistryType = rtype;
	}

	/**
	 * Should we use the "standard" set of handlers for the visualization?
	 * If false, use dummy (stub) handlers.
	 * \see IOptions.h
	 */
	bool useVizStandardHandlers() const override {
		return vizStandardHandlersFlag;
	}

	/**
	 * Set the vizStandardHandlersFlag.
	 * \see IOptions.h
	 */
	void setVizStandardHandlers(bool flag) override {
		vizStandardHandlersFlag = flag;
	}

	/**
	 * Obtain the name of the material to be used for simulation.
	 * \see IOptions.h
	 */
	std::string getMaterial() const override {
		return materialName;
	}

	/**
	 * Set the name of the material to be used for the simulation.
	 * \see IOptions.h
	 */
	void setMaterial(const std::string& material) override {
		materialName = material;
	}

	/**
	 * Obtain the value of the concentration for the vacancies.
	 * \see IOptions.h
	 */
	double getInitialVConcentration() const override {
		return initialVConcentration;
	}

	/**
	 * Set the value of the concentration for the vacancies.
	 * \see IOptions.h
	 */
	void setInitialVConcentration(double conc) override {
		initialVConcentration = conc;
	}

	/**
	 * Obtain the number of dimensions for the simulation.
	 * \see IOptions.h
	 */
	int getDimensionNumber() const override {
		return dimensionNumber;
	}

	/**
	 * Set the number of dimensions for the simulation.
	 * \see IOptions.h
	 */
	void setDimensionNumber(int number) override {
		dimensionNumber = number;
	}

	/**
	 * Obtain the value of the void portion for the simulation.
	 * \see IOptions.h
	 */
	double getVoidPortion() const override {
		return voidPortion;
	}

	/**
	 * Set the value of the void portion for the surface to grow.
	 * \see IOptions.h
	 */
	void setVoidPortion(double portion) override {
		voidPortion = portion;
	}

	/**
	 * Should we use a regular grid on the x direction?
	 * \see IOptions.h
	 */
	bool useRegularXGrid() const override {
		return useRegularGridFlag;
	}

	/**
	 * Set the useRegularGridFlag.
	 * \see IOptions.h
	 */
	void setRegularXGrid(bool flag) override {
		useRegularGridFlag = flag;
	}

	/**
	 * Obtain the physical process map.
	 *
	 * @return The map
	 */
	std::map<std::string, bool> getProcesses() const override {
		return processMap;
	}

	/**
	 * Set the physical process map.
	 *
	 * @param map The map
	 */
	void setProcesses(std::map<std::string, bool> map) override {
		processMap = map;
	}

	/**
	 * Obtain the string listing the wanted GB.
	 * \see IOptions.h
	 */
	std::string getGbString() const override {
		return gbList;
	}

	/**
	 * Set the string listing the wanted GB.
	 * \see IOptions.h
	 */
	void setGbString(const std::string& gbString) override {
		gbList = gbString;
	}

	/**
	 * Obtain the minimum size for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingMin() const override {
		return groupingMin;
	}

	/**
	 * Set the minimum size for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingMin(int size) override {
		groupingMin = size;
	}

	/**
	 * Obtain the first width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthA() const override {
		return groupingWidthA;
	}

	/**
	 * Set the first width for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingWidthA(int width) override {
		groupingWidthA = width;
	}

	/**
	 * Obtain the second width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthB() const override {
		return groupingWidthB;
	}

	/**
	 * Set the second width for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingWidthB(int width) override {
		groupingWidthB = width;
	}

	/**
	 * Obtain the value of the intensity of the sputtering yield to be used.
	 * \see IOptions.h
	 */
	double getSputteringYield() const override {
		return sputteringYield;
	}

	/**
	 * Set the value for the sputtering yield to use.
	 * \see IOptions.h
	 */
	void setSputteringYield(double yield) override {
		sputteringYield = yield;
	}

	/**
	 * To know if we should use the HDF5 file.
	 * \see IOptions.h
	 */
	bool useHDF5() const override {
		return useHDF5Flag;
	}

	/**
	 * Set the useHDF5Flag.
	 * \see IOptions.h
	 */
	void setHDF5Flag(bool flag) override {
		useHDF5Flag = flag;
	}

	/**
	 * To know if we should use the phase cut.
	 * \see IOptions.h
	 */
	bool usePhaseCut() const override {
		return usePhaseCutFlag;
	}

	/**
	 * Set the usePhaseCutFlag.
	 * \see IOptions.h
	 */
	void setPhaseCutFlag(bool flag) override {
		usePhaseCutFlag = flag;
	}

	/**
	 * Obtain the maximum value of impurities (He or Xe) to be used.
	 * \see IOptions.h
	 */
	int getMaxImpurity() const override {
		return maxImpurity;
	}

	/**
	 * Set the maximum value of impurities to use.
	 * \see IOptions.h
	 */
	void setMaxImpurity(int max) override {
		maxImpurity = max;
	}

	/**
	 * Obtain the maximum value of vacancies to be used.
	 * \see IOptions.h
	 */
	int getMaxV() const override {
		return maxV;
	}

	/**
	 * Set the maximum value of vacancies to use.
	 * \see IOptions.h
	 */
	void setMaxV(int max) override {
		maxV = max;
	}

	/**
	 * Obtain the maximum value of interstitials to be used.
	 * \see IOptions.h
	 */
	int getMaxI() const override {
		return maxI;
	}

	/**
	 * Set the maximum value of interstitials to use.
	 * \see IOptions.h
	 */
	void setMaxI(int max) override {
		maxI = max;
	}

	/**
	 * Obtain the number of grid points in the depth direction to be used.
	 * \see IOptions.h
	 */
	int getNX() const override {
		return nX;
	}

	/**
	 * Set the number of grid points in the depth direction to use.
	 * \see IOptions.h
	 */
	void setNX(int n) override {
		nX = n;
	}

	/**
	 * Obtain the value of the step size in the depth direction to be used.
	 * \see IOptions.h
	 */
	double getXStepSize() const override {
		return xStepSize;
	}

	/**
	 * Set the value for the step size in the depth direction to use.
	 * \see IOptions.h
	 */
	void setXStepSize(double stepSize) override {
		xStepSize = stepSize;
	}

	/**
	 * Obtain the number of grid points in the Y direction to be used.
	 * \see IOptions.h
	 */
	int getNY() const override {
		return nY;
	}

	/**
	 * Set the number of grid points in the Y direction to use.
	 * \see IOptions.h
	 */
	void setNY(int n) override {
		nY = n;
	}

	/**
	 * Obtain the value of the step size in the Y direction to be used.
	 * \see IOptions.h
	 */
	double getYStepSize() const override {
		return yStepSize;
	}

	/**
	 * Set the value for the step size in the Y direction to use.
	 * \see IOptions.h
	 */
	void setYStepSize(double stepSize) override {
		yStepSize = stepSize;
	}

	/**
	 * Obtain the number of grid points in the Z direction to be used.
	 * \see IOptions.h
	 */
	int getNZ() const override {
		return nZ;
	}

	/**
	 * Set the number of grid points in the Z direction to use.
	 * \see IOptions.h
	 */
	void setNZ(int n) override {
		nZ = n;
	}

	/**
	 * Obtain the value of the step size in the Z direction to be used.
	 * \see IOptions.h
	 */
	double getZStepSize() const override {
		return zStepSize;
	}

	/**
	 * Set the value for the step size in the Z direction to use.
	 * \see IOptions.h
	 */
	void setZStepSize(double stepSize) override {
		zStepSize = stepSize;
	}

	/**
	 * Indicate whether to output the reaction network to a file
     * (e.g., to support debugging the network).
	 *
     * @param shouldWriteNetwork Whether to write reaction network 
     *        to file once created.
	 * @param fileName Name of the file to which the network should be written.
     *        If fileName = "-", network is written to standard output.
     *        Ignored unless shouldWriteNetwork is true.
	 */
	virtual void setNetworkDebugOptions(bool _shouldWriteNetwork = false,
                                        std::string fileName = "network.txt") override {
        shouldWriteDebugNetwork = _shouldWriteNetwork;
        debugNetworkFilename = fileName;
    }


	/**
	 * Retrieve user's settings for network debugging.
	 *
     * @return Pair (b, f) where b indicates whether to write the 
     *          reaction network to a file, and f indicates the filename to use.
	 */
    virtual std::pair<bool, std::string> getNetworkDebugOptions() const override {
        return std::make_pair(shouldWriteDebugNetwork, debugNetworkFilename);
    }
};
//end class Options

} /* namespace xolotlCore */

#endif // OPTIONS_H
