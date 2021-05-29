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

import java.io.IOException;
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
    private static int STOP_LOOPS = 1;

    /** The time limit in seconds required for termination. */
    private static int STOP_TIME = 7200;

    /**
     * Creates the directories for the result files.
     */
    public static void createDirectories() {
	try {
	    Files.createDirectories(Paths.get("Results/APSA/"));
	    Files.createDirectories(Paths.get("Results/AVA/"));
	    Files.createDirectories(Paths.get("Results/HAPSA/"));
	    Files.createDirectories(Paths.get("Results/APSA/"));
	    Files.createDirectories(Paths.get("Results/HAPSA/"));
	} catch (IOException e) {
	    System.err.println("Directory could not be created in Main.createDirectories().");
	}
    }

    /**
     * Runs the specified methods.
     * 
     * @param args no arguments required
     */
    public static void main(String[] args) {
	System.out.println("Simulating a store with 50 shelves and 400 products...");
	Store store = (new StoreSimulator(400, 50, 0)).simulate();

	long startTime = System.currentTimeMillis() / 1000;
	Solution apsa = solveAPSA(store);
	long endTime = System.currentTimeMillis() / 1000;
	apsa.writeSolution("Results/APSA/50_400.txt", endTime - startTime);

	double[] gammas = new double[] { 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10 };
	for (double gamma : gammas) {
	    startTime = System.currentTimeMillis() / 1000;
	    Solution hapsa = solveHAPSA(store, gamma, 0);
	    endTime = System.currentTimeMillis() / 1000;
	    hapsa.writeSolution("Results/HAPSA/50_400_" + gamma + "_0.txt", endTime - startTime);
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
		ParameterSet params = initModel.readParameterSet("Parameters/SSP/" + obj + "_" + ni + "_" + nj);
		initModel.setParameterSet(params);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (w[pr] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && w[nonSelected.indexOf(j1)] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (nonSelected.contains(j2) && w[nonSelected.indexOf(j2)] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
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
		ParameterSet params = initModel.readParameterSet("Parameters/SSP/APSA_" + ni + "_" + nj);
		initModel.setParameterSet(params);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (w[pr] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && w[nonSelected.indexOf(j1)] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (nonSelected.contains(j2) && w[nonSelected.indexOf(j2)] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
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
		ParameterSet params = initModel.readParameterSet("Parameters/SSP/HAPSA_" + ni + "_" + nj);
		initModel.setParameterSet(params);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int pr = 0; pr < w.length; pr++) {
		    int j = nonSelected.get(pr);
		    if (w[pr] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[pr];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (nonSelected.contains(j1) && w[nonSelected.indexOf(j1)] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (nonSelected.contains(j2) && w[nonSelected.indexOf(j2)] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
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
	    ParameterSet params = continuous.readParameterSet("Parameters/CONT/" + obj + "_" + ni + "_" + nj);
	    continuous.setParameterSet(params);
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
		return Double.compare(obj1, obj2);
	    }
	}

	int loops = 0;
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
		    ParameterSet params = model.readParameterSet("Parameters/HAPSA/" + obj + "_" + ni + "_" + nj);
		    model.setParameterSet(params);
		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.002);
		    double maxTime = 0.000015 * Math.pow((ni + nj), 2.8082);
		    model.setParam(IloCplex.Param.TimeLimit, maxTime);
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
		    if (newObjective <= objective) {
			loops += N_REOPT / store.getShelves().size();
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
	    ParameterSet params = continuous.readParameterSet("Parameters/CONT/APSA_" + ni + "_" + nj);
	    continuous.setParameterSet(params);
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
		return Double.compare(obj1, obj2);
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
		Store partialStore = incumbent.partialStore(selected, consideredProducts);
		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    model.setObjective(Model.Objective.APSA);

		    ParameterSet params = model.readParameterSet("Parameters/HAPSA/APSA_" + ni + "_" + nj);
		    model.setParameterSet(params);
		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.002);
		    double maxTime = 0.000015 * Math.pow((ni + nj), 2.8082);
		    model.setParam(IloCplex.Param.TimeLimit, maxTime);
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
		    if (newObjective <= objective) {
			loops += N_REOPT / store.getShelves().size();
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
	    ParameterSet params = continuous.readParameterSet("Parameters/CONT/HAPSA_" + ni + "_" + nj);
	    continuous.setParameterSet(params);
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
		return Double.compare(obj1, obj2);
	    }
	}

	int loops = 0;
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
		Store partialStore = incumbent.partialStore(selected, consideredProducts);
		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    model.setGamma(gamma);
		    model.setTheta(theta);
		    model.setObjective(Model.Objective.HAPSA);
		    ParameterSet params = model.readParameterSet("Parameters/HAPSA/HAPSA_" + ni + "_" + nj);
		    model.setParameterSet(params);
		    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.002);
		    double maxTime = 0.000015 * Math.pow((ni + nj), 2.8082);
		    model.setParam(IloCplex.Param.TimeLimit, maxTime);
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
		    if (newObjective <= objective) {
			loops += N_REOPT / store.getShelves().size();
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
}
