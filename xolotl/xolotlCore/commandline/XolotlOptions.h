#ifndef XOLOTLOPTIONS_H
#define XOLOTLOPTIONS_H

#include <string>
#include "Options.h"

namespace xolotlCore {

class XolotlOptions : public Options
{
private:

    // Use the "W" (tungsten) set of handlers?
    bool useWHandlers;

    // Use the constant temperature set of handlers?
    bool useConstTempHandlers;

    // Use the temperature profile set of handlers?
    bool useTempProfileHandlers;

    // Use the "standard" set of handlers?
    bool useStdHandlers;

    // Name of the input network file.
    std::string netFileName;

    // Value of the constant temperature in Kelvin
    double constTemp;

    // Name of the input temperature profile file.
    std::string tempProfileFileName;

    // Deal with the material handler selection option.
    // @param arg Argument given to the material handler selection option.
    bool handleMaterialOption( std::string arg );

    // Callback when have seen the material handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleMaterialOptionCB( Options* opts, std::string arg );

    // Deal with the constant temperature handler selection option.
    // @param arg Argument given to the constant temperature handler selection option.
    bool handleConstTemperatureOption( std::string arg);

    // Callback when have seen the constant temperature handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleConstTemperatureOptionCB( Options* opts, std::string arg);

    // Deal with the temperature file handler selection option.
    // @param arg Argument given to the temperature file handler selection option.
    bool handleTemperatureFileOption( std::string arg);

    // Callback when have seen the temperature file handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleTemperatureFileOptionCB( Options* opts, std::string arg);

    // Deal with the handler selection option.
    // @param arg Argument given to the handler selection option.
    bool handleHandlersOption( std::string arg );

    // Callback when have seen the handlers option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handleHandlersOptionCB( Options* opts, std::string arg );


    // Deal with the PETSc argument delimiter option.
    // @param arg Unused
    bool handlePetscOption( std::string arg );

    // Callback when have seen the PETSc argument delimiter option.
    // @param opts Options object for the program.
    // @param arg Argument provided to the option.
    static bool handlePetscOptionCB( Options* opts, std::string arg );

public:
    XolotlOptions( void );


    // Parse the given command line for user-configurable settings.
    // We assume that the executable file/path has been skipped before
    // calling this method.  (E.g., the program's main() function
    // called this with something like 
    //   xopts.parseCommandLine( argc - 1, argv + 1 );
    //
    // @param argc The number of arguments in the argv vector.
    // @param argv Vector of argument strings.
    // @return Number of command line arguments used.
    virtual int parseCommandLine( int argc, char* argv[] );


    // Show our help message.
    // @param os The output stream upon which to print the help message.
    virtual void showHelp( std::ostream& os ) const;

    // Should we use "W" (tungsten) handlers?
    // If false, use Fe (iron) handlers.
    // TODO will we ever have more than {W,Fe} such that
    // we will need an enum?
    // @return true if program should use tungsten handlers, false if
    // should use iron handlers.
    bool useTungstenHandlers( void ) const  { return useWHandlers; }

    // Should we use const temperature handlers?
    // @return true if program should use const temp handlers
    bool useConstTemperatureHandlers( void ) const  { return useConstTempHandlers; }

    // Should we use temperature profile handlers?
    // @return true if program should use temperature profile handlers
    bool useTemperatureProfileHandlers( void ) const  { return useTempProfileHandlers; }

    // Obtain the name of the file containing the temperature profile data.
    // @return Name of the temperature file.
    std::string getTempProfileFilename( void ) const    { return tempProfileFileName; }

    // Obtain the value of the constant temperature to be used.
    // @return Constant temperature.
    double getConstTemperature( void ) const  { return constTemp; }

    // Should we use the "standard" set of handlers?
    // If false, use dummy (stub) handlers.
    // TODO will we ever have more than {dummy, standard} such that 
    // we will need an enum?
    // @return true if program should use standard handlers, false if 
    // should use dummy handlers.
    bool useStandardHandlers( void ) const  { return useStdHandlers; }


    // Obtain the name of the file holding the input network.
    // @return Name of the input network file.
    std::string getNetworkFilename( void ) const    { return netFileName; }
};

};

#endif // XOLOTLOPTIONS_H
