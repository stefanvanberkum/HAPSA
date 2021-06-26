/*
 * Main.java
 * 
 * v1.0
 *
 * 14 May 2021
 *
 * Copyright (c) 2021 Stefan van Berkum
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
 * OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package hapsa;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.concurrent.ThreadLocalRandom;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.ParameterSet;

/**
 * The Main class provides the main execution environment, and contains methods
 * for the optimization-based heuristic approach.
 *
 * @author Stefan van Berkum
 *
 */
public class Main {

    /**
     * The required optimality gap for each iteration of the re-optimization
     * procedure.
     **/
    private static double MIP_GAP = 0.001;

    /** The number of shelves to be selected for each re-optimization run. */
    private static int N_REOPT = 4;

    /**
     * The relative gap between the upper bound and the incumbent solution required
     * for termination.
     */
    private static double STOP_GAP = 0.005;

    /**
     * The number of loops over all shelves without improvement to the incumbent
     * solution required for termination.
     */
    private static int STOP_LOOPS = 10;

    /** The time limit in seconds required for termination. */
    private static int STOP_TIME = 18000;

    /** The time limit on each iteration of the re-optimization procedure. **/
    private static int TIME_LIMIT = 120;

    /**
     * True if tuned parameters should be used. If set to false, default CPLEX
     * parameters will be used.
     **/
    private static boolean USE_PARAM = true;

    /**
     * Creates the directories for the result files.
     */
    public static void createDirectories() {
	try {
	    String[] types = new String[] { "APSA", "AVA", "HAPSA", "HLUR", "VIS", "APSA_2D" };
	    for (String type : types) {
		Files.createDirectories(Paths.get("Results/" + type + "/Summary/"));
		Files.createDirectories(Paths.get("Results/" + type + "/Variables/"));
		Files.createDirectories(Paths.get("Results/" + type + "/CSV/"));
	    }
	    Files.createDirectories(Paths.get("Log/"));
	} catch (IOException e) {
	    System.err.println("Directory could not be created in Main.createDirectories().");
	}
    }

    /**
     * Runs the specified methods.
     * 
     * NOTE. Be sure to save the main_log file whenever you start a new run (the
     * file will be overwritten).
     * 
     * @param args no arguments required
     */
    public static void main(String[] args) {
	createDirectories();
	PrintStream console = System.out;
	try {
	    PrintStream log = new PrintStream(new File("Log/main_log.txt"));
	    System.setOut(log);
	} catch (FileNotFoundException e) {
	    System.err.println("File not found in Main.main(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	// Solve a 3D implementation of APSA (as in our paper).
	System.out.println("Simulating a store with 30 shelves and 240 products...");
	Store store = (new StoreSimulator(240, 30, 0)).simulate();

	console.println("Running APSA...");
	long startTime = System.currentTimeMillis() / 1000;
	Solution apsa = solveAPSA(store);
	long endTime = System.currentTimeMillis() / 1000;
	apsa.writeSolution("Results/APSA/", "30_240", endTime - startTime);
	console.println("APSA ran for " + (endTime - startTime) + " seconds.");

	double[] gammas = new double[] { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
		105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200 };
	for (double gamma : gammas) {
	    console.println("Running APSA with visibility penalty (gamma = " + gamma + ")...");
	    startTime = System.currentTimeMillis() / 1000;
	    Solution hapsa = solve(store, Model.Objective.VIS, gamma);
	    endTime = System.currentTimeMillis() / 1000;
	    hapsa.writeSolution("Results/VIS/", "30_240_" + gamma, endTime - startTime);
	    console.println("APSA with visibility penalty (gamma = " + gamma + ") ran for " + (endTime - startTime)
		    + " seconds.");
	}

	double gamma = 100;
	double[] thetas = new double[] { 0.0001, 0.00025, 0.0005, 0.001 };
	for (double theta : thetas) {
	    console.println("Running HAPSA (gamma = " + gamma + ", theta = " + theta + ")...");
	    startTime = System.currentTimeMillis() / 1000;
	    Solution hapsa = solveHAPSA(store, gamma, theta);
	    endTime = System.currentTimeMillis() / 1000;
	    hapsa.writeSolution("Results/HAPSA/", "30_240_" + gamma + "_" + theta, endTime - startTime);
	    console.println("HAPSA (gamma = " + gamma + ", theta = " + theta + ") ran for " + (endTime - startTime)
		    + " seconds.");
	}

	// Solve a 2D implementation of APSA (as in original paper).
	// Note that changing the static variables in this way is not recommended,
	// it is just for the sake of demonstration.
	double CAPACITY = 6;
	double[] END_SEG = new double[] { 0.06, 0.1 };
	double HEALTH_SCORE_LB = 1;
	double HEALTH_SCORE_UB = 100;
	int HORIZONTAL = 3;
	double[] HORIZONTAL_CAT = new double[] { 0.05, 0.25, 0.45, 0.65, 0.85 };
	double MAX_PROFIT_LB = 1;
	double MAX_PROFIT_UB = 25;
	double MAX_SPACE_UB = 6;
	double[] MIDDLE_SEG = new double[] { 0.0, 0.05 };
	double MIN_ALLOCATED = 0.1;
	double MIN_SPACE_LB = 1;
	double MIN_SPACE_UB = 3;
	double L_FRAC = 0.0;
	double H1_FRAC = 0.0;
	double H2_FRAC = 0.0;
	double H3_FRAC = 0.0;
	int VERTICAL = 1;
	double[] VERTICAL_CAT = new double[] { 1.0 };
	StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL, HORIZONTAL_CAT,
		MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED, MIN_SPACE_LB, MIN_SPACE_UB,
		L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);

	int runs = 3;
	int[] shelves = new int[] { 50, 100 };
	int[] products = new int[] { 400, 800 };
	int[] ts = new int[] { 2, 3, 4 };
	double[] epsilons = new double[] { 0.015, 0.01, 0.005 };
	STOP_TIME = 3600;

	for (int i = 0; i < shelves.length; i++) {
	    for (int j = 0; j < runs; j++) {
		L_FRAC = 0.0;
		H1_FRAC = 0.0;
		H2_FRAC = 0.0;
		H3_FRAC = 0.0;
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		store = (new StoreSimulator(products[i], shelves[i], j)).simulate();

		N_REOPT = 4;
		STOP_GAP = 0.005;

		// Regular APSA.
		console.println("Running APSA (shelves = " + shelves[i] + ", products = " + products[i] + ")...");

		startTime = System.currentTimeMillis() / 1000;
		try (HAPSA model = new HAPSA(store)) {
		    model.setObjective(Model.Objective.APSA);
		    ParameterSet params = model.readParameterSet("Parameters/CONT/APSA_30_240");
		    model.setParameterSet(params);
		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, STOP_GAP);
		    model.setParam(IloCplex.Param.TimeLimit, STOP_TIME);
		    model.solve();

		    double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
		    double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
		    double[][] y = new double[store.getSegments().size()][store.getProducts().size()];

		    for (int sh = 0; sh < store.getShelves().size(); sh++) {
			x[sh] = model.getValues(model.getX()[sh]);
		    }

		    for (int seg = 0; seg < store.getSegments().size(); seg++) {
			s[seg] = model.getValues(model.getS()[seg]);
			y[seg] = model.getValues(model.getY()[seg]);
		    }

		    Solution result = new Solution(s, x, y, store);
		    result.updateObjectiveAPSA();
		    double upperBound = model.getObjValue() * (1 + model.getMIPRelativeGap())
			    + Math.pow(10, -10) * model.getMIPRelativeGap();
		    result.setUpperBound(upperBound);
		    endTime = System.currentTimeMillis() / 1000;
		    result.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_" + j,
			    endTime - startTime);
		    console.println("APSA (shelves = " + shelves[i] + ", products = " + products[i] + ") ran for "
			    + (endTime - startTime) + " seconds.");
		} catch (IloException e) {
		    System.out.println("APSA could not be solved in Main.main(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}

		// Heuristic with affinities.
		// With five sets from L.
		L_FRAC = (double) 5 / products[i];
		H1_FRAC = 0.0;
		H2_FRAC = 0.0;
		H3_FRAC = 0.0;
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		Store lStore = (new StoreSimulator(products[i], shelves[i], j)).simulate();
		console.println(
			"Running heuristic with L (shelves = " + shelves[i] + ", products = " + products[i] + ")...");
		startTime = System.currentTimeMillis() / 1000;
		Solution model = solveAPSA(lStore);
		endTime = System.currentTimeMillis() / 1000;
		model.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_L_" + j,
			endTime - startTime);
		console.println("Heuristic with L (shelves = " + shelves[i] + ", products = " + products[i]
			+ ") ran for " + (endTime - startTime) + " seconds.");

		// With five sets from H1.
		L_FRAC = 0.0;
		H1_FRAC = (double) 5 / products[i];
		H2_FRAC = 0.0;
		H3_FRAC = 0.0;
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		Store h1Store = (new StoreSimulator(products[i], shelves[i], j)).simulate();
		console.println(
			"Running heuristic with H1 (shelves = " + shelves[i] + ", products = " + products[i] + ")...");
		startTime = System.currentTimeMillis() / 1000;
		model = solveAPSA(h1Store);
		endTime = System.currentTimeMillis() / 1000;
		model.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_H1_" + j,
			endTime - startTime);
		console.println("Heuristic with H1 (shelves = " + shelves[i] + ", products = " + products[i]
			+ ") ran for " + (endTime - startTime) + " seconds.");

		// With five sets from H2.
		L_FRAC = 0.0;
		H1_FRAC = 0.0;
		H2_FRAC = (double) 5 / products[i];
		H3_FRAC = 0.0;
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		Store h2Store = (new StoreSimulator(products[i], shelves[i], j)).simulate();
		console.println(
			"Running heuristic with H2 (shelves = " + shelves[i] + ", products = " + products[i] + ")...");
		startTime = System.currentTimeMillis() / 1000;
		model = solveAPSA(h2Store);
		endTime = System.currentTimeMillis() / 1000;
		model.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_H2_" + j,
			endTime - startTime);
		console.println("Heuristic with H2 (shelves = " + shelves[i] + ", products = " + products[i]
			+ ") ran for " + (endTime - startTime) + " seconds.");

		// With five sets from H3.
		L_FRAC = 0.0;
		H1_FRAC = 0.0;
		H2_FRAC = 0.0;
		H3_FRAC = (double) 5 / products[i];
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		Store h3Store = (new StoreSimulator(products[i], shelves[i], j)).simulate();
		console.println(
			"Running heuristic with H3 (shelves = " + shelves[i] + ", products = " + products[i] + ")...");
		startTime = System.currentTimeMillis() / 1000;
		model = solveAPSA(h3Store);
		endTime = System.currentTimeMillis() / 1000;
		model.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_H3_" + j,
			endTime - startTime);
		console.println("Heuristic with H3 (shelves = " + shelves[i] + ", products = " + products[i]
			+ ") ran for " + (endTime - startTime) + " seconds.");

		// With five sets from L, H1, H2, and H3.
		L_FRAC = (double) 5 / products[i];
		H1_FRAC = (double) 5 / products[i];
		H2_FRAC = (double) 5 / products[i];
		H3_FRAC = (double) 5 / products[i];
		StoreSimulator.setParameters(CAPACITY, END_SEG, HEALTH_SCORE_LB, HEALTH_SCORE_UB, HORIZONTAL,
			HORIZONTAL_CAT, MAX_PROFIT_LB, MAX_PROFIT_UB, MAX_SPACE_UB, MIDDLE_SEG, MIN_ALLOCATED,
			MIN_SPACE_LB, MIN_SPACE_UB, L_FRAC, H1_FRAC, H2_FRAC, H3_FRAC, VERTICAL, VERTICAL_CAT);
		Store affinityStore = (new StoreSimulator(products[i], shelves[i], j)).simulate();
		console.println("Running heuristic with affinities (shelves = " + shelves[i] + ", products = "
			+ products[i] + ")...");
		startTime = System.currentTimeMillis() / 1000;
		model = solveAPSA(affinityStore);
		endTime = System.currentTimeMillis() / 1000;
		model.writeSolution("Results/APSA_2D/", shelves[i] + "_" + products[i] + "_affinities_" + j,
			endTime - startTime);
		console.println("Heuristic with affinities (shelves = " + shelves[i] + ", products = " + products[i]
			+ ") ran for " + (endTime - startTime) + " seconds.");

		// Heuristic without affinities.
		for (int t : ts) {
		    for (double epsilon : epsilons) {
			N_REOPT = t;
			STOP_GAP = epsilon;

			console.println("Running heuristic (shelves = " + shelves[i] + ", products = " + products[i]
				+ ", t = " + t + ", epsilon = " + epsilon + ")...");
			startTime = System.currentTimeMillis() / 1000;
			model = solveAPSA(store);
			endTime = System.currentTimeMillis() / 1000;
			model.writeSolution("Results/APSA_2D/",
				shelves[i] + "_" + products[i] + "_" + t + "_" + epsilon + "_" + j,
				endTime - startTime);
			console.println("Heuristic (shelves = " + shelves[i] + ", products = " + products[i] + ", t = "
				+ t + ", epsilon = " + epsilon + ") ran for " + (endTime - startTime) + " seconds.");
		    }
		}
	    }
	}
    }

    /**
     * Solves the model using the optimization-based heuristic approach for
     * objective functions of type: AVA, HLUR, or VIS.
     * 
     * @param store the store to be considered
     * @param obj   the objective function type, one of: AVA, HLUR, or VIS
     * @param param the parameter that corresponds to the objective function type
     * @return the solution
     */
    public static Solution solve(Store store, Model.Objective obj, double param) {
	Solution solution = initialize(store, obj, param);
	return reOptimize(solution, store, obj, param);
    }

    /**
     * Solves the model using the optimization-based heuristic approach for the
     * objective function type APSA.
     * 
     * @param store the store to be considered
     * @return the solution
     */
    public static Solution solveAPSA(Store store) {
	Solution solution = initializeAPSA(store);
	return reOptimizeAPSA(solution, store);
    }

    /**
     * Solves the model using the optimization-based heuristic approach for the
     * objective function type HAPSA.
     * 
     * @param store the store to be considered
     * @param gamma the gamma parameter for the visibility penalty
     * @param theta the theta parameter for the healthy-left, unhealthy-right
     *              approach
     * @return the solution
     */
    public static Solution solveHAPSA(Store store, double gamma, double theta) {
	Solution solution = initializeHAPSA(store, gamma, theta);
	return reOptimizeHAPSA(solution, store, gamma, theta);
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for objective functions of type: AVA, HLUR, or VIS.
     * 
     * @param store the store to be considered
     * @param obj   the objective function type, one of: AVA, HLUR, or VIS
     * @param param the parameter that corresponds to the objective function type
     * @return the solution
     * @throws IllegalStateException when the algorithm cannot find a feasible
     *                               solution for a shelf
     */
    private static Solution initialize(Store store, Model.Objective obj, double param) throws IllegalStateException {
	ArrayList<Shelf> sortedShelves = new ArrayList<Shelf>(store.getShelves());
	Collections.sort(sortedShelves);
	HashSet<Product> selected = new HashSet<Product>();
	double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
	double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
	double[][] y = new double[store.getSegments().size()][store.getProducts().size()];
	Solution result = new Solution(s, x, y, store);

	System.out.println("Initiating " + obj + " initialization procedure...");

	// Loop through all shelves.
	for (int sh = 0; sh < sortedShelves.size(); sh++) {
	    int i = store.getShelves().indexOf(sortedShelves.get(sh));
	    Shelf shelf = sortedShelves.get(sh);
	    ArrayList<Shelf> shelves = new ArrayList<Shelf>(1);
	    shelves.add(shelf);

	    TreeSet<Integer> products = new TreeSet<Integer>();
	    for (Product product : store.getProducts()) {
		if (!selected.contains(product)) {
		    products.add(store.getProducts().indexOf(product));
		}
	    }
	    ArrayList<Integer> nonSelected = new ArrayList<Integer>(products);

	    Store oneShelfStore = result.partialStore(shelves, nonSelected);
	    try (SSP initModel = new SSP(oneShelfStore)) {
		System.out.println("Shelf " + (sh + 1) + " out of " + sortedShelves.size() + ".");

		switch (obj) {
		case AVA:
		    initModel.setLambda(param);
		    break;
		case HLUR:
		    initModel.setTheta(param);
		    break;
		default:
		    initModel.setGamma(param);
		    break;
		}
		initModel.setObjective(obj);
		int ni = store.getShelves().size();
		int nj = store.getProducts().size();

		if (USE_PARAM) {
		    ParameterSet params = initModel.readParameterSet("Parameters/SSP/" + obj + "_" + ni + "_" + nj);
		    initModel.setParameterSet(params);
		}

		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (Math.round(w[pr]) == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && Math.round(w[nonSelected.indexOf(j1)]) == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
			updateRelationships(store, selected, j2);
		    }
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int seg = 0; seg < oneShelfStore.getSegments().size(); seg++) {
		    Segment segment = oneShelfStore.getSegments().get(seg);
		    int k = store.getSegments().indexOf(segment);

		    for (int pr = 0; pr < oneShelfStore.getProducts().size(); pr++) {
			Product product = oneShelfStore.getProducts().get(pr);
			int j = store.getProducts().indexOf(product);

			s[k][j] = initModel.getValue(sInit[seg][pr]);
			y[k][j] = initModel.getValue(yInit[seg][pr]);
		    }
		}

		if (selected.size() == store.getProducts().size()) {
		    break;
		}
	    } catch (IloException e) {
		System.err.println("SSP model could not be created or solved in Main.initialize(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for the objective function of type APSA.
     * 
     * @param store the store to be considered
     * @return the solution
     * @throws IllegalStateException when the algorithm cannot find a feasible
     *                               solution for a shelf
     */
    private static Solution initializeAPSA(Store store) throws IllegalStateException {
	ArrayList<Shelf> sortedShelves = new ArrayList<Shelf>(store.getShelves());
	Collections.sort(sortedShelves);
	HashSet<Product> selected = new HashSet<Product>();
	double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
	double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
	double[][] y = new double[store.getSegments().size()][store.getProducts().size()];
	Solution result = new Solution(s, x, y, store);

	System.out.println("Initiating APSA initialization procedure...");

	// Loop through all shelves.
	for (int sh = 0; sh < sortedShelves.size(); sh++) {
	    int i = store.getShelves().indexOf(sortedShelves.get(sh));
	    Shelf shelf = sortedShelves.get(sh);
	    ArrayList<Shelf> shelves = new ArrayList<Shelf>(1);
	    shelves.add(shelf);

	    TreeSet<Integer> products = new TreeSet<Integer>();
	    for (Product product : store.getProducts()) {
		if (!selected.contains(product)) {
		    products.add(store.getProducts().indexOf(product));
		}
	    }
	    ArrayList<Integer> nonSelected = new ArrayList<Integer>(products);

	    Store oneShelfStore = result.partialStore(shelves, nonSelected);
	    try (SSP initModel = new SSP(oneShelfStore)) {
		System.out.println("Shelf " + (sh + 1) + " out of " + sortedShelves.size() + ".");

		initModel.setObjective(Model.Objective.APSA);
		int ni = store.getShelves().size();
		int nj = store.getProducts().size();

		if (USE_PARAM) {
		    ParameterSet params = initModel.readParameterSet("Parameters/SSP/APSA_" + ni + "_" + nj);
		    initModel.setParameterSet(params);
		}

		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (Math.round(w[pr]) == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && Math.round(w[nonSelected.indexOf(j1)]) == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
			updateRelationships(store, selected, j2);
		    }
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int seg = 0; seg < oneShelfStore.getSegments().size(); seg++) {
		    Segment segment = oneShelfStore.getSegments().get(seg);
		    int k = store.getSegments().indexOf(segment);

		    for (int pr = 0; pr < oneShelfStore.getProducts().size(); pr++) {
			Product product = oneShelfStore.getProducts().get(pr);
			int j = store.getProducts().indexOf(product);

			s[k][j] = initModel.getValue(sInit[seg][pr]);
			y[k][j] = initModel.getValue(yInit[seg][pr]);
		    }
		}

		if (selected.size() == store.getProducts().size()) {
		    break;
		}
	    } catch (IloException e) {
		System.err.println("SSP model could not be created or solved in Main.initializeAPSA(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for the objective function of type HAPSA.
     * 
     * @param store the store to be considered
     * @param gamma the gamma parameter for the visibility penalty
     * @param theta the theta parameter for the healthy-left, unhealthy-right
     *              approach
     * @return the solution
     * @throws IllegalStateException when the algorithm cannot find a feasible
     *                               solution for a shelf
     */
    private static Solution initializeHAPSA(Store store, double gamma, double theta) throws IllegalStateException {
	ArrayList<Shelf> sortedShelves = new ArrayList<Shelf>(store.getShelves());
	Collections.sort(sortedShelves);
	HashSet<Product> selected = new HashSet<Product>();
	double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
	double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
	double[][] y = new double[store.getSegments().size()][store.getProducts().size()];
	Solution result = new Solution(s, x, y, store);

	System.out.println("Initiating HAPSA initialization procedure...");

	// Loop through all shelves.
	for (int sh = 0; sh < sortedShelves.size(); sh++) {
	    int i = store.getShelves().indexOf(sortedShelves.get(sh));
	    Shelf shelf = sortedShelves.get(sh);
	    ArrayList<Shelf> shelves = new ArrayList<Shelf>(1);
	    shelves.add(shelf);

	    TreeSet<Integer> products = new TreeSet<Integer>();
	    for (Product product : store.getProducts()) {
		if (!selected.contains(product)) {
		    products.add(store.getProducts().indexOf(product));
		}
	    }
	    ArrayList<Integer> nonSelected = new ArrayList<Integer>(products);

	    Store oneShelfStore = result.partialStore(shelves, nonSelected);
	    try (SSP initModel = new SSP(oneShelfStore)) {
		System.out.println("Shelf " + (sh + 1) + " out of " + sortedShelves.size() + ".");

		initModel.setGamma(gamma);
		initModel.setTheta(theta);
		initModel.setObjective(Model.Objective.HAPSA);
		int ni = store.getShelves().size();
		int nj = store.getProducts().size();

		if (USE_PARAM) {
		    ParameterSet params = initModel.readParameterSet("Parameters/SSP/HAPSA_" + ni + "_" + nj);
		    initModel.setParameterSet(params);
		}

		initModel.setParam(IloCplex.Param.TimeLimit, TIME_LIMIT);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (Math.round(w[pr]) == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && Math.round(w[nonSelected.indexOf(j1)]) == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
			updateRelationships(store, selected, j2);
		    }
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j2) && Math.round(w[nonSelected.indexOf(j2)]) == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
			updateRelationships(store, selected, j1);
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int seg = 0; seg < oneShelfStore.getSegments().size(); seg++) {
		    Segment segment = oneShelfStore.getSegments().get(seg);
		    int k = store.getSegments().indexOf(segment);

		    for (int pr = 0; pr < oneShelfStore.getProducts().size(); pr++) {
			Product product = oneShelfStore.getProducts().get(pr);
			int j = store.getProducts().indexOf(product);

			s[k][j] = initModel.getValue(sInit[seg][pr]);
			y[k][j] = initModel.getValue(yInit[seg][pr]);
		    }
		}

		if (selected.size() == store.getProducts().size()) {
		    break;
		}
	    } catch (IloException e) {
		System.err.println("SSP model could not be created or solved in Main.initializeHAPSA(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Re-optimizes the shelves iteratively by means of the MIP-based
     * re-optimization procedure, for objective function types: AVA, HLUR, or VIS.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     * @param obj   the objective function type, one of: AVA, HLUR, or VIS
     * @param param the parameter that corresponds to the objective function type
     * @return the solution
     */
    private static Solution reOptimize(Solution init, Store store, Model.Objective obj, double param) {
	System.out.println("Initiating " + obj + " re-optimization procedure...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Solve continuous relaxation of the model.
	double upperBound = Double.MAX_VALUE;
	try (HAPSA continuous = new HAPSA(store)) {
	    switch (obj) {
	    case AVA:
		continuous.setLambda(param);
		break;
	    case HLUR:
		continuous.setTheta(param);
		break;
	    default:
		continuous.setGamma(param);
		break;
	    }
	    continuous.setObjective(obj);
	    continuous.relax();

	    if (USE_PARAM) {
		ParameterSet params = continuous.readParameterSet("Parameters/CONT/" + obj + "_" + ni + "_" + nj);
		continuous.setParameterSet(params);
	    }

	    boolean feasible = continuous.solve();
	    if (!feasible) {
		throw new IllegalStateException("No feasible upper bound found.");
	    }
	    upperBound = continuous.getObjValue();
	    System.out.println("Upper bound: " + upperBound + " (" + continuous.getStatus() + ")");
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in Main.reOptimize(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
	incumbent.setUpperBound(upperBound);
	double objective = incumbent.updateObjective(obj, param);
	System.out.println("Objective = " + objective);

	/**
	 * The ObjectiveComparator class compares shelves based on their current
	 * objective value contribution.
	 *
	 * @author Stefan van Berkum
	 *
	 */
	class ObjectiveComparator implements Comparator<Shelf> {
	    @Override
	    public int compare(Shelf shelf1, Shelf shelf2) {
		double obj1 = incumbent.getShelfObjective(shelf1, obj, param);
		double obj2 = incumbent.getShelfObjective(shelf2, obj, param);
		return Double.compare(obj2, obj1);
	    }
	}

	double loops = 0;
	long startTime = System.currentTimeMillis() / 1000;
	long time = System.currentTimeMillis() / 1000;

	while (((upperBound - objective) / upperBound > STOP_GAP) && (loops < STOP_LOOPS)
		&& (time - startTime < STOP_TIME)) {
	    TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	    shelves.addAll(store.getShelves());

	    while (shelves.size() > store.getShelves().size() % N_REOPT) {
		int nPerLevel = shelves.size() / N_REOPT;
		ArrayList<Shelf> selected = new ArrayList<Shelf>(N_REOPT);
		Iterator<Shelf> iter = shelves.iterator();

		// Select a shelf from each level.
		for (int k = 1; k < N_REOPT + 1; k++) {
		    ArrayList<Shelf> level;
		    int randIndex;
		    if (k < N_REOPT) {
			level = new ArrayList<Shelf>(nPerLevel);
			for (int i = 0; i < nPerLevel; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nPerLevel);
		    } else {
			// Dump remaining shelves into last level.
			int nLast = shelves.size() - (N_REOPT - 1) * nPerLevel;
			level = new ArrayList<Shelf>(nLast);
			for (int i = 0; i < nLast; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nLast);
		    }
		    selected.add(level.get(randIndex));
		}
		shelves.removeAll(selected);

		ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);

		// Remove all products j from (j, j') that have asymmetric assortment affinity,
		// when j' is selected on a shelf that is not considered.
		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		}

		// Remove all products j (j') from (j, j') that have allocation affinity, when
		// j' (j) is selected on a shelf that is not considered.
		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		    if (consideredProducts.contains(j2) && !consideredProducts.contains(j1)) {
			// Product 2 cannot be considered in this run as product 1 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j2));
			updateRelationships(store, consideredProducts, j2);
		    }
		}

		Store partialStore = incumbent.partialStore(selected, consideredProducts);

		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    switch (obj) {
		    case AVA:
			model.setLambda(param);
			break;
		    case HLUR:
			model.setTheta(param);
			break;
		    default:
			model.setGamma(param);
			break;
		    }
		    model.setObjective(obj);

		    if (USE_PARAM) {
			ParameterSet params = model.readParameterSet("Parameters/HAPSA/" + obj + "_" + ni + "_" + nj);
			model.setParameterSet(params);
		    }

		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
		    model.setParam(IloCplex.Param.TimeLimit, TIME_LIMIT);
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    // Collect partial results and update incumbent solution.
		    int[] shelfIndices = new int[N_REOPT];
		    int nseg = partialStore.getSegments().size();
		    int npr = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][npr];
		    double[][] partialS = new double[nseg][npr];
		    double[][] partialY = new double[nseg][npr];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nseg; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double newObjective = incumbent.updateObjective(obj, param);
		    if (newObjective - objective < 0.01) {
			loops += (double) N_REOPT / store.getShelves().size();
		    } else {
			loops = 0;
		    }
		    objective = newObjective;
		    System.out.println("Objective = " + objective);
		} catch (IloException e) {
		    System.err.println("Subproblem could not be solved in Main.reOptimize(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	    time = System.currentTimeMillis() / 1000;
	    System.out.println("Elapsed time re-optimization: " + (time - startTime) + " seconds");
	}
	if ((upperBound - objective) / upperBound <= STOP_GAP) {
	    incumbent.setTerminationReason(Solution.Termination.GAP);
	} else if (loops >= STOP_LOOPS) {
	    incumbent.setTerminationReason(Solution.Termination.LOOP);
	} else {
	    incumbent.setTerminationReason(Solution.Termination.TIME);
	}
	return incumbent;
    }

    /**
     * Re-optimizes the shelves iteratively by means of the MIP-based
     * re-optimization procedure, for the objective function type APSA.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     * @return the solution
     */
    private static Solution reOptimizeAPSA(Solution init, Store store) {
	System.out.println("Initiating APSA re-optimization procedure...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Solve continuous relaxation of the model.
	double upperBound = Double.MAX_VALUE;
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setObjective(Model.Objective.APSA);
	    continuous.relax();

	    if (USE_PARAM) {
		ParameterSet params = continuous.readParameterSet("Parameters/CONT/APSA_" + ni + "_" + nj);
		continuous.setParameterSet(params);
	    }

	    boolean feasible = continuous.solve();
	    if (!feasible) {
		throw new IllegalStateException("No feasible upper bound found.");
	    }
	    upperBound = continuous.getObjValue();
	    System.out.println("Upper bound: " + upperBound + " (" + continuous.getStatus() + ")");
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in Main.reOptimize(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
	incumbent.setUpperBound(upperBound);
	double objective = incumbent.updateObjectiveAPSA();
	System.out.println("Objective = " + objective);

	/**
	 * The ObjectiveComparator class compares shelves based on their current
	 * objective value contribution.
	 *
	 * @author Stefan van Berkum
	 *
	 */
	class ObjectiveComparator implements Comparator<Shelf> {
	    @Override
	    public int compare(Shelf shelf1, Shelf shelf2) {
		double obj1 = incumbent.getShelfObjectiveAPSA(shelf1);
		double obj2 = incumbent.getShelfObjectiveAPSA(shelf2);
		return Double.compare(obj2, obj1);
	    }
	}

	double loops = 0;
	long startTime = System.currentTimeMillis() / 1000;
	long time = System.currentTimeMillis() / 1000;

	while (((upperBound - objective) / upperBound > STOP_GAP) && (loops < STOP_LOOPS)
		&& (time - startTime < STOP_TIME)) {

	    TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	    shelves.addAll(store.getShelves());

	    while (shelves.size() > store.getShelves().size() % N_REOPT) {
		int nPerLevel = shelves.size() / N_REOPT;
		ArrayList<Shelf> selected = new ArrayList<Shelf>(N_REOPT);
		Iterator<Shelf> iter = shelves.iterator();

		// Select a shelf from each level.
		for (int k = 1; k < N_REOPT + 1; k++) {
		    ArrayList<Shelf> level;
		    int randIndex;
		    if (k < N_REOPT) {
			level = new ArrayList<Shelf>(nPerLevel);
			for (int i = 0; i < nPerLevel; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nPerLevel);
		    } else {
			// Dump remaining shelves into last level.
			int nLast = shelves.size() - (N_REOPT - 1) * nPerLevel;
			level = new ArrayList<Shelf>(nLast);
			for (int i = 0; i < nLast; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nLast);
		    }
		    selected.add(level.get(randIndex));
		}
		shelves.removeAll(selected);

		ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);

		// Remove all products j from (j, j') that have asymmetric assortment affinity,
		// when j' is selected on a shelf that is not considered.
		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		}

		// Remove all products j (j') from (j, j') that have allocation affinity, when
		// j' (j) is selected on a shelf that is not considered.
		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		    if (consideredProducts.contains(j2) && !consideredProducts.contains(j1)) {
			// Product 2 cannot be considered in this run as product 1 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j2));
			updateRelationships(store, consideredProducts, j2);
		    }
		}

		Store partialStore = incumbent.partialStore(selected, consideredProducts);

		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    model.setObjective(Model.Objective.APSA);

		    if (USE_PARAM) {
			ParameterSet params = model.readParameterSet("Parameters/HAPSA/APSA_" + ni + "_" + nj);
			model.setParameterSet(params);
		    }

		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
		    model.setParam(IloCplex.Param.TimeLimit, TIME_LIMIT);
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    // Collect partial results and update incumbent solution.
		    int[] shelfIndices = new int[N_REOPT];
		    int nseg = partialStore.getSegments().size();
		    int npr = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][npr];
		    double[][] partialS = new double[nseg][npr];
		    double[][] partialY = new double[nseg][npr];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nseg; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double newObjective = incumbent.updateObjectiveAPSA();

		    if (newObjective - objective < 0.01) {
			loops += (double) N_REOPT / store.getShelves().size();
		    } else {
			loops = 0;
		    }
		    objective = newObjective;
		    System.out.println("Objective = " + objective);
		} catch (IloException e) {
		    System.err.println("Subproblem could not be solved in Main.reOptimize(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	    time = System.currentTimeMillis() / 1000;
	    System.out.println("Elapsed time re-optimization: " + (time - startTime) + " seconds");
	}
	if ((upperBound - objective) / upperBound <= STOP_GAP) {
	    incumbent.setTerminationReason(Solution.Termination.GAP);
	} else if (loops >= STOP_LOOPS) {
	    incumbent.setTerminationReason(Solution.Termination.LOOP);
	} else {
	    incumbent.setTerminationReason(Solution.Termination.TIME);
	}
	return incumbent;
    }

    /**
     * Re-optimizes the shelves iteratively by means of the MIP-based
     * re-optimization procedure, for the objective function type APSA.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     * @param gamma the gamma parameter for the visibility penalty
     * @param theta the theta parameter for the healthy-left, unhealthy-right
     *              approach
     * @return the solution
     */
    private static Solution reOptimizeHAPSA(Solution init, Store store, double gamma, double theta) {
	System.out.println("Initiating HAPSA re-optimization procedure...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Solve continuous relaxation of the model.
	double upperBound = Double.MAX_VALUE;
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setGamma(gamma);
	    continuous.setTheta(theta);
	    continuous.setObjective(Model.Objective.HAPSA);
	    continuous.relax();

	    if (USE_PARAM) {
		ParameterSet params = continuous.readParameterSet("Parameters/CONT/HAPSA_" + ni + "_" + nj);
		continuous.setParameterSet(params);
	    }

	    boolean feasible = continuous.solve();
	    if (!feasible) {
		throw new IllegalStateException("No feasible upper bound found.");
	    }
	    upperBound = continuous.getObjValue();
	    System.out.println("Upper bound: " + upperBound + " (" + continuous.getStatus() + ")");
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in Main.reOptimize(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
	incumbent.setUpperBound(upperBound);
	double objective = incumbent.updateObjectiveHAPSA(gamma, theta);
	System.out.println("Objective = " + objective);

	/**
	 * The ObjectiveComparator class compares shelves based on their current
	 * objective value contribution.
	 *
	 * @author Stefan van Berkum
	 *
	 */
	class ObjectiveComparator implements Comparator<Shelf> {
	    @Override
	    public int compare(Shelf shelf1, Shelf shelf2) {
		double obj1 = incumbent.getShelfObjectiveHAPSA(shelf1, gamma, theta);
		double obj2 = incumbent.getShelfObjectiveHAPSA(shelf2, gamma, theta);
		return Double.compare(obj2, obj1);
	    }
	}

	double loops = 0;
	long startTime = System.currentTimeMillis() / 1000;
	long time = System.currentTimeMillis() / 1000;

	while (((upperBound - objective) / upperBound > STOP_GAP) && (loops < STOP_LOOPS)
		&& (time - startTime < STOP_TIME)) {
	    TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	    shelves.addAll(store.getShelves());

	    while (shelves.size() > store.getShelves().size() % N_REOPT) {
		int nPerLevel = shelves.size() / N_REOPT;
		ArrayList<Shelf> selected = new ArrayList<Shelf>(N_REOPT);
		Iterator<Shelf> iter = shelves.iterator();

		// Select a shelf from each level.
		for (int k = 1; k < N_REOPT + 1; k++) {
		    ArrayList<Shelf> level;
		    int randIndex;
		    if (k < N_REOPT) {
			level = new ArrayList<Shelf>(nPerLevel);
			for (int i = 0; i < nPerLevel; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nPerLevel);
		    } else {
			// Dump remaining shelves into last level.
			int nLast = shelves.size() - (N_REOPT - 1) * nPerLevel;
			level = new ArrayList<Shelf>(nLast);
			for (int i = 0; i < nLast; i++) {
			    level.add(iter.next());
			}
			randIndex = ThreadLocalRandom.current().nextInt(0, nLast);
		    }
		    selected.add(level.get(randIndex));
		}
		shelves.removeAll(selected);

		ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);

		// Remove all products j from (j, j') that have asymmetric assortment affinity,
		// when j' is selected on a shelf that is not considered.
		for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		}

		// Remove all products j (j') from (j, j') that have allocation affinity, when
		// j' (j) is selected on a shelf that is not considered.
		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (consideredProducts.contains(j1) && !consideredProducts.contains(j2)) {
			// Product 1 cannot be considered in this run as product 2 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j1));
			updateRelationships(store, consideredProducts, j1);
		    }
		    if (consideredProducts.contains(j2) && !consideredProducts.contains(j1)) {
			// Product 2 cannot be considered in this run as product 1 is already placed on
			// a shelf that is not considered.
			consideredProducts.remove(consideredProducts.indexOf(j2));
			updateRelationships(store, consideredProducts, j2);
		    }
		}

		Store partialStore = incumbent.partialStore(selected, consideredProducts);

		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    model.setGamma(gamma);
		    model.setTheta(theta);
		    model.setObjective(Model.Objective.HAPSA);

		    if (USE_PARAM) {
			ParameterSet params = model.readParameterSet("Parameters/HAPSA/HAPSA_" + ni + "_" + nj);
			model.setParameterSet(params);
		    }

		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
		    model.setParam(IloCplex.Param.TimeLimit, TIME_LIMIT);
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    // Collect partial results and update incumbent solution.
		    int[] shelfIndices = new int[N_REOPT];
		    int nseg = partialStore.getSegments().size();
		    int npr = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][npr];
		    double[][] partialS = new double[nseg][npr];
		    double[][] partialY = new double[nseg][npr];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nseg; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double newObjective = incumbent.updateObjectiveHAPSA(gamma, theta);
		    if (newObjective - objective < 0.01) {
			loops += (double) N_REOPT / store.getShelves().size();
		    } else {
			loops = 0;
		    }
		    objective = newObjective;
		    System.out.println("Objective = " + objective);
		} catch (IloException e) {
		    System.err.println("Subproblem could not be solved in Main.reOptimize(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	    time = System.currentTimeMillis() / 1000;
	    System.out.println("Elapsed time re-optimization: " + (time - startTime) + " seconds");
	}
	if ((upperBound - objective) / upperBound <= STOP_GAP) {
	    incumbent.setTerminationReason(Solution.Termination.GAP);
	} else if (loops >= STOP_LOOPS) {
	    incumbent.setTerminationReason(Solution.Termination.LOOP);
	} else {
	    incumbent.setTerminationReason(Solution.Termination.TIME);
	}
	return incumbent;
    }

    /**
     * Updates the affinity relationships for the re-optimization algorithm when a
     * product is removed from consideration. For a product j that is removed from
     * consideration, it removes any product j' from consideration as well if (j,
     * j') or (j', j) in H1 (symmetric assortment affinity) or if (j', j) in H2
     * (asymmetric assortment affinity).
     * 
     * @param store              the store
     * @param consideredProducts the list of currently considered products
     * @param removedIndex       the index of the product that is removed from
     *                           consideration
     */
    private static void updateRelationships(Store store, ArrayList<Integer> consideredProducts, int removedIndex) {
	// For a product j that is removed from consideration, remove all products j'
	// from consideration if (j, j') or (j', j) in H1.
	for (SymmetricPair pair : store.getSymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (consideredProducts.contains(j1) && removedIndex == j2) {
		// Product 1 cannot be considered in this run as product 2 is also not
		// considered.
		consideredProducts.remove(consideredProducts.indexOf(j1));

		// Recursively update other relations.
		updateRelationships(store, consideredProducts, j1);
	    }
	    if (consideredProducts.contains(j2) && removedIndex == j1) {
		// Product 2 cannot be considered in this run as product 1 is also not
		// considered.
		consideredProducts.remove(consideredProducts.indexOf(j2));

		// Recursively update other relations.
		updateRelationships(store, consideredProducts, j2);
	    }
	}

	// For a product j that is removed from consideration, remove all products j'
	// from consideration if (j', j) in H2.
	for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (consideredProducts.contains(j1) && removedIndex == j2) {
		// Product 1 cannot be considered in this run as product 2 is also not
		// considered.
		consideredProducts.remove(consideredProducts.indexOf(j1));

		// Recursively update other relations.
		updateRelationships(store, consideredProducts, j1);
	    }
	}
    }

    /**
     * Updates the affinity relationships for the initialization algorithm when a
     * product is removed from consideration. For a product j that is removed from
     * consideration, it removes any product j' from consideration as well if (j,
     * j') or (j', j) in H1 (symmetric assortment affinity) or if (j', j) in H2
     * (asymmetric assortment affinity).
     * 
     * @param store        the store
     * @param selected     the list of currently selected products (which are not
     *                     further considered)
     * @param removedIndex the index of the product that is removed from
     *                     consideration
     */
    private static void updateRelationships(Store store, HashSet<Product> selected, int removedIndex) {
	// For a product j that is removed from consideration, remove all products j'
	// from consideration if (j, j') or (j', j) in H1.
	for (SymmetricPair pair : store.getSymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (!selected.contains(store.getProducts().get(j1)) && removedIndex == j2) {
		// Product 1 cannot be considered in subsequent runs as product 2 is also not
		// considered.
		selected.add(store.getProducts().get(j1));

		// Recursively update other relations.
		updateRelationships(store, selected, j1);
	    }
	    if (!selected.contains(store.getProducts().get(j2)) && removedIndex == j1) {
		// Product 2 cannot be considered in subsequent runs as product 1 is also not
		// considered.
		selected.add(store.getProducts().get(j2));

		// Recursively update other relations.
		updateRelationships(store, selected, j2);
	    }
	}

	// For a product j that is removed from consideration, remove all products j'
	// from consideration if (j', j) in H2.
	for (AsymmetricPair pair : store.getAsymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (!selected.contains(store.getProducts().get(j1)) && removedIndex == j2) {
		// Product 1 cannot be considered in subsequent runs as product 2 is also not
		// considered.
		selected.add(store.getProducts().get(j1));

		// Recursively update other relations.
		updateRelationships(store, selected, j1);
	    }
	}
    }
}
