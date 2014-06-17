package gov.ornl.xolotl.preprocessor.test;

import static org.junit.Assert.*;

import org.junit.Test;

import gov.ornl.xolotl.preprocessor.Arguments;

import uk.co.flamingpenguin.jewel.cli.ArgumentValidationException;
import uk.co.flamingpenguin.jewel.cli.CliFactory;

/**
 * This class is responsible for testing the Preprocessor class
 */
public class ArgumentsTest {

	/**
	 * This operation checks the default Arguments.
	 */
	@Test
	public void testDefaultArguments() throws ArgumentValidationException {

		// Local Declarations
		final Arguments args = CliFactory.parseArguments(Arguments.class, new String [] {});

		// Check that the default material is W
		assertEquals("W", args.getMaterial());

		// Check that the default startTemp is 1000
		assertEquals("1000", args.getStartTemp());

		// Check if there is a default tempFile 
		assertEquals(false, args.isTempFile());

		// Check that the default heFlux is 2.5e27
		assertEquals("2.5e27", args.getHeFlux());
		
		// Check that the default perfHandler is dummy
		assertEquals("dummy", args.getPerfHandler());
		
		// Check that the default vizHandler is dummy
		assertEquals("dummy", args.getVizHandler());
		
		// Check that the default checkpoint is true
		assertEquals("true", args.getCheckpoint());
		
		// Check that the default networkFile is test.h5
		assertEquals("test.h5", args.getNetworkFile());
		
		// Check the default petscArgs
		assertEquals("-da_grid_x 10 -ts_final_time 1000"
				+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200"
				+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant"
				+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
				args.getPetscArgs());
	}
	
	/**
	 * This operation tests that default parameter values are only overridden
	 * if they are specified via the command line
	 */
	@Test
	public void testSpecifiedArguments() throws ArgumentValidationException {

		// Local Declarations
		final Arguments args = CliFactory.parseArguments(Arguments.class, 
				new String [] {"--startTemp", "900", "--material", "Fe",
								"--perfHandler", "std"});

		// Check that the default material is W
		assertEquals("Fe", args.getMaterial());

		// Check that the default startTemp is 1000
		assertEquals("900", args.getStartTemp());

		// Check if there is a default tempFile 
		assertEquals(false, args.isTempFile());

		// Check that the default heFlux is 2.5e27
		assertEquals("2.5e27", args.getHeFlux());
		
		// Check that the default perfHandler is dummy
		assertEquals("std", args.getPerfHandler());
		
		// Check that the default vizHandler is dummy
		assertEquals("dummy", args.getVizHandler());
		
		// Check that the default checkpoint is true
		assertEquals("true", args.getCheckpoint());
		
		// Check that the default networkFile is test.h5
		assertEquals("test.h5", args.getNetworkFile());
		
		// Check the default petscArgs
		assertEquals("-da_grid_x 10 -ts_final_time 1000"
				+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200"
				+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant"
				+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor",
				args.getPetscArgs());
	}

}
