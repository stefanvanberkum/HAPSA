/*
 * ParameterTuner.java
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
 * The ParamterTuner class provides methods to tune the CPLEX parameters for the
 * optimization-based heuristic approach to shelf space planning. The optimal
 * parameters are determined for the first iteration of the initialization
 * procedure, the continuous relaxation, and the first iteration of the
 * re-optimization procedure. The results are saved to a file and automatically
 * loaded into the {@link Main} class.
 *
 * NOTE. Re-running the parameter tuning overwrites the previous parameter
 * tuning files.
 *
 * @author Stefan van Berkum
 *
 */
public class ParameterTuner {

    /**
     * The required optimality gap for each iteration of the re-optimization
     * procedure.
     **/
    private static double MIP_GAP = 0.001;

    /** The number of shelves to be selected for each re-optimization run. */
    private static int N_REOPT = 4;

    /**
     * Creates the directories for the parameter optimization result files.
     */
    public static void createDirectories() {
	try {
	    Files.createDirectories(Paths.get("Parameters/SSP/"));
	    Files.createDirectories(Paths.get("Parameters/CONT/"));
	    Files.createDirectories(Paths.get("Parameters/HAPSA/"));
	} catch (IOException e) {
	    System.err.println("Directory could not be created in ParamterTuner.createDirectories().");
	}
    }

    /**
     * Runs the parameter tuning.
     * 
     * @param args no arguments required
     */
    public static void main(String[] args) {
	System.out.println("Creating directories...");
	createDirectories();
	System.out.println("Simulating a store with 30 shelves and 240 products...");
	Store store = (new StoreSimulator(240, 30, 0)).simulate();
	System.out.println("Tuning APSA...");
	tuneAPSA(store);
	System.out.println("Tuning APSA with visibility penalty...");
	tune(store, Model.Objective.VIS, 100);
	System.out.println("Tuning HAPSA...");
	tuneHAPSA(store, 100, 5);
    }

    /**
     * Tunes the CPLEX parameters for objective functions of type: AVA, HLUR, or
     * VIS.
     * 
     * @param store the store to be considered
     * @param obj   the objective function type, one of: AVA, HLUR, or VIS
     * @param param the parameter that corresponds to the objective function type
     */
    public static void tune(Store store, Model.Objective obj, double param) {
	Solution solution = initialize(store, obj, param);
	tuneReOpt(solution, store, obj, param);
    }

    /**
     * Tunes the CPLEX parameters for the objective function type APSA.
     * 
     * @param store the store to be considered
     */
    public static void tuneAPSA(Store store) {
	Solution solution = initializeAPSA(store);
	tuneReOptAPSA(solution, store);
    }

    /**
     * Tunes the CPLEX parameters for the objective function type HAPSA.
     * 
     * @param store the store to be considered
     * @param gamma the gamma parameter for the visibility penalty
     * @param theta the theta parameter for the healthy-left, unhealthy-right
     *              approach
     */
    public static void tuneHAPSA(Store store, double gamma, double theta) {
	Solution solution = initializeHAPSA(store, gamma, theta);
	tuneReOptHAPSA(solution, store, gamma, theta);
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for objective functions of type: AVA, HLUR, or VIS. Additionally, this
     * method tunes the CPLEX parameters for the first iteration and saves the
     * results to a file.
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

		if (sh == 0) {
		    // Collect tuning information for the first shelf.
		    initModel.tuneParam();
		    ParameterSet params = initModel.getParameterSet();
		    initModel.writeParameterSet(params, "Parameters/SSP/" + obj + "_" + ni + "_" + nj);
		} else {
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
		System.err.println("SSP model could not be created or solved in ParameterTuner.initialize(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for the objective function of type APSA. Additionally, this method
     * tunes the CPLEX parameters for the first iteration and saves the results to a
     * file.
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

		if (sh == 0) {
		    // Collect tuning information for the first shelf.
		    initModel.tuneParam();
		    ParameterSet params = initModel.getParameterSet();
		    initModel.writeParameterSet(params, "Parameters/SSP/APSA_" + ni + "_" + nj);
		} else {
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
		System.err.println("SSP model could not be created or solved in ParameterTuner.initializeAPSA(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Initializes the model by optimizing each shelf separately using the SSP(i*)
     * model, for the objective function of type HAPSA. Additionally, this method
     * tunes the CPLEX parameters for the first iteration and saves the results to a
     * file.
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

		if (sh == 0) {
		    // Collect tuning information for the first shelf.
		    initModel.tuneParam();
		    ParameterSet params = initModel.getParameterSet();
		    initModel.writeParameterSet(params, "Parameters/SSP/HAPSA_" + ni + "_" + nj);
		} else {
		    ParameterSet params = initModel.readParameterSet("Parameters/SSP/HAPSA_" + ni + "_" + nj);
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
		System.err.println("SSP model could not be created or solved in ParameterTuner.initializeHAPSA(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
	return result;
    }

    /**
     * Tunes the CPLEX parameters for the MIP-based re-optimization procedure, for
     * objective function types: AVA, HLUR, or VIS. The parameters are tuned for the
     * continuous relaxation and for the first iteration of the re-optimization
     * procedure, and the results are saved to a file.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     * @param obj   the objective function type, one of: AVA, HLUR, or VIS
     * @param param the parameter that corresponds to the objective function type
     */
    private static void tuneReOpt(Solution init, Store store, Model.Objective obj, double param) {
	System.out.println("Initiating " + obj + " re-optimization tuning...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Tune continuous relaxation parameters of the model.
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

	    // Collect tuning information.
	    continuous.tuneParam();
	    ParameterSet params = continuous.getParameterSet();
	    continuous.writeParameterSet(params, "Parameters/CONT/" + obj + "_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in ParameterTuner.tuneReOpt(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
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

	TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	shelves.addAll(store.getShelves());

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

	    // Collect tuning information for the re-optimization run.
	    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
	    model.tuneParam();
	    ParameterSet params = model.getParameterSet();
	    model.writeParameterSet(params, "Parameters/HAPSA/" + obj + "_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.err.println("Subproblem could not be solved in ParameterTuner.tuneReOpt(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Tunes the CPLEX parameters for the MIP-based re-optimization procedure, for
     * the objective function type APSA. The parameters are tuned for the continuous
     * relaxation and for the first iteration of the re-optimization procedure, and
     * the results are saved to a file.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     */
    private static void tuneReOptAPSA(Solution init, Store store) {
	System.out.println("Initiating APSA re-optimization tuning...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Solve continuous relaxation of the model.
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setObjective(Model.Objective.APSA);
	    continuous.relax();

	    // Collect tuning information.
	    continuous.tuneParam();
	    ParameterSet params = continuous.getParameterSet();
	    continuous.writeParameterSet(params, "Parameters/CONT/APSA_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in ParameterTuner.tuneReOptAPSA(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
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

	TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	shelves.addAll(store.getShelves());

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

	ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);
	Store partialStore = incumbent.partialStore(selected, consideredProducts);
	try (HAPSA model = new HAPSA(partialStore)) {
	    model.initializePartial(incumbent, consideredProducts);
	    model.setObjective(Model.Objective.APSA);

	    // Collect tuning information for the re-optimization run.
	    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
	    model.tuneParam();
	    ParameterSet params = model.getParameterSet();
	    model.writeParameterSet(params, "Parameters/HAPSA/APSA_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.err.println("Subproblem could not be solved in ParameterTuner.tuneReOptAPSA(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Tunes the CPLEX parameters for the MIP-based re-optimization procedure, for
     * the objective function type APSA. The parameters are tuned for the continuous
     * relaxation and for the first iteration of the re-optimization procedure, and
     * the results are saved to a file.
     * 
     * @param init  the initial solution provided by the initialization procedure
     * @param store the store to be considered
     * @param gamma the gamma parameter for the visibility penalty
     * @param theta the theta parameter for the healthy-left, unhealthy-right
     *              approach
     */
    private static void tuneReOptHAPSA(Solution init, Store store, double gamma, double theta) {
	System.out.println("Initiating HAPSA re-optimization tuning...");
	int ni = store.getShelves().size();
	int nj = store.getProducts().size();

	// Solve continuous relaxation of the model.
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setGamma(gamma);
	    continuous.setTheta(theta);
	    continuous.setObjective(Model.Objective.HAPSA);
	    continuous.relax();

	    // Collect tuning information.
	    continuous.tuneParam();
	    ParameterSet params = continuous.getParameterSet();
	    continuous.writeParameterSet(params, "Parameters/CONT/HAPSA_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in ParameterTuner.tuneReOptHAPSA(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
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

	TreeSet<Shelf> shelves = new TreeSet<Shelf>(new ObjectiveComparator());
	shelves.addAll(store.getShelves());

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

	ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);
	Store partialStore = incumbent.partialStore(selected, consideredProducts);
	try (HAPSA model = new HAPSA(partialStore)) {
	    model.initializePartial(incumbent, consideredProducts);
	    model.setGamma(gamma);
	    model.setTheta(theta);
	    model.setObjective(Model.Objective.HAPSA);

	    // Collect tuning information for the re-optimization run.
	    model.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, MIP_GAP);
	    model.tuneParam();
	    ParameterSet params = model.getParameterSet();
	    model.writeParameterSet(params, "Parameters/HAPSA/HAPSA_" + ni + "_" + nj);
	} catch (IloException e) {
	    System.err.println("Subproblem could not be solved in ParameterTuner.tuneReOptHAPSA(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }
}
