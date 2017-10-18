#ifndef NETWORKDEBUGOPTIONHANDLER_H
#define NETWORKDEBUGOPTIONHANDLER_H

// Includes
#include <string>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * NetworkDebugOptionHandler handles the option to generate
 * the network if not loaded from the HDF5 file.
 */
class NetworkDebugOptionHandler : public OptionHandler {
public:

    /**
     * The default constructor
     */
    NetworkDebugOptionHandler() :
        OptionHandler("netDebug",
                    "netDebug                          "
                    "Whether to output debugging information for reaction networks.\n"
                    "Format: shouldWriteNetwork [fileName]\n"
                    "    where shouldWriteNetwork = [true | false] (default false)\n"
                    "    and fileName is name of file to write to\n"
                    "        (default \"network.txt\", use \"-\" for standard output"
                ) {
    }

    /**
     * The destructor
     */
    ~NetworkDebugOptionHandler() {
    }

    /**
     * This method will set the IOptions network parameters
     * to the values given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The list of PETSc options.
     */
    bool handler(IOptions *opt, const std::string& arg) {
        // Set the flag to not use the HDF5 file
        opt->setHDF5Flag(false);

        // Build an input stream from the argument string.
        xolotlCore::TokenizedLineReader<std::string> reader;
        auto argSS = std::make_shared<std::istringstream>(arg);
        reader.setInputStream(argSS);

        // Break the argument into tokens.
        auto tokens = reader.loadLine();

        // Interpret the tokens.
        std::string networkFilename("network.txt");
        bool shouldWriteNetworkFile = (tokens[0] == "true");

        if(shouldWriteNetworkFile and tokens.size() > 1) {
            networkFilename = tokens[1];            
        }
        opt->setNetworkDebugOptions(shouldWriteNetworkFile, networkFilename);

        return true;
    }

};

} /* namespace xolotlCore */

#endif // NETWORKDEBUGOPTIONHANDLER_H
