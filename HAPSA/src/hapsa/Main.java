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
    private static int STOP_TIME = 3600;

    public static void main(String[] args) {
	System.out.println("Simulating the first store...");
	Store store1 = (new StoreSimulator(240, 30, 1)).simulate();
	Solution apsa1 = solveAPSA(store1);

	// Store store2 = (new StoreSimulator(320, 40, 2)).simulate();
	// Store store3 = (new StoreSimulator(400, 50, 3)).simulate();
	// Store store4 = (new StoreSimulator(480, 60, 4)).simulate();
	// Store store5 = (new StoreSimulator(800, 100, 5)).simulate();
    }

    public static Solution solve(Store store, SSP.Objective obj, double param) {
	Solution solution = initialize(store, obj, param);
	switch (obj) {
	case AVA:
	    return reOptimize(solution, store, HAPSA.Objective.AVA, param);
	case HLUR:
	    return reOptimize(solution, store, HAPSA.Objective.HLUR, param);
	default:
	    return reOptimize(solution, store, HAPSA.Objective.VIS, param);
	}
    }

    public static Solution solveAPSA(Store store) {
	Solution solution = initializeAPSA(store);
	return reOptimizeAPSA(solution, store);
    }

    public static Solution solveHAPSA(Store store, double gamma, double theta) {
	Solution solution = initializeHAPSA(store, gamma, theta);
	return reOptimizeHAPSA(solution, store, gamma, theta);
    }

    private static Solution initialize(Store store, SSP.Objective obj, double param) throws IllegalStateException {
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
		System.out.println("Shelf " + sh + " out of " + sortedShelves.size() + ".");

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
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int j = 0; j < w.length; j++) {
		    if (w[j] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[j];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (w[j1] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (w[j2] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int k = 0; k < sInit.length; k++) {
		    s[i * sInit.length + k] = initModel.getValues(sInit[k]);
		    y[i * sInit.length + k] = initModel.getValues(yInit[k]);
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
	return new Solution(s, x, y, store);
    }

    private static Solution initializeAPSA(Store store) throws IllegalStateException {
	ArrayList<Shelf> sortedShelves = new ArrayList<Shelf>(store.getShelves());
	Collections.sort(sortedShelves);
	HashSet<Product> selected = new HashSet<Product>();
	double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
	double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
	double[][] y = new double[store.getSegments().size()][store.getProducts().size()];

	System.out.println("Initiating APSA initialization procedure...");

	for (int sh = 0; sh < sortedShelves.size(); sh++) {
	    int i = store.getShelves().indexOf(sortedShelves.get(sh));
	    try (SSP initModel = new SSP(store, sortedShelves.get(sh), selected)) {
		System.out.println("Shelf " + sh + " out of " + sortedShelves.size() + ".");

		initModel.setObjective(SSP.Objective.APSA);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int j = 0; j < w.length; j++) {
		    if (w[j] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[j];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (w[j1] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (w[j2] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int k = 0; k < sInit.length; k++) {
		    s[i * sInit.length + k] = initModel.getValues(sInit[k]);
		    y[i * sInit.length + k] = initModel.getValues(yInit[k]);
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
	return new Solution(s, x, y, store);
    }

    private static Solution initializeHAPSA(Store store, double gamma, double theta) throws IllegalStateException {
	ArrayList<Shelf> sortedShelves = new ArrayList<Shelf>(store.getShelves());
	Collections.sort(sortedShelves);
	HashSet<Product> selected = new HashSet<Product>();
	double[][] s = new double[store.getSegments().size()][store.getProducts().size()];
	double[][] x = new double[store.getShelves().size()][store.getProducts().size()];
	double[][] y = new double[store.getSegments().size()][store.getProducts().size()];

	System.out.println("Initiating HAPSA initialization procedure...");

	for (int sh = 0; sh < sortedShelves.size(); sh++) {
	    int i = store.getShelves().indexOf(sortedShelves.get(sh));
	    try (SSP initModel = new SSP(store, sortedShelves.get(sh), selected)) {
		System.out.println("Shelf " + sh + " out of " + sortedShelves.size() + ".");

		initModel.setGamma(gamma);
		initModel.setTheta(theta);
		initModel.setObjective(SSP.Objective.HAPSA);
		boolean feasible = initModel.solve();

		if (!feasible) {
		    throw new IllegalStateException("No feasible solution found.");
		}

		System.out.println("Feasible solution found! Type: " + initModel.getStatus() + "\n");

		double[] w = initModel.getValues(initModel.getW());
		for (int j = 0; j < w.length; j++) {
		    if (w[j] == 1) {
			// This product is selected on this shelf.
			selected.add(store.getProducts().get(j));
		    }
		    x[i][j] = w[j];
		}

		for (SymmetricPair pair : store.getAllocationAffinity()) {
		    int j1 = pair.getIndex1();
		    int j2 = pair.getIndex2();
		    if (w[j1] == 1) {
			// Product 2 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j2));
		    }
		    if (w[j2] == 1) {
			// Product 1 cannot be placed on another shelf.
			selected.add(store.getProducts().get(j1));
		    }
		}

		IloNumVar[][] sInit = initModel.getS();
		IloIntVar[][] yInit = initModel.getY();
		for (int k = 0; k < sInit.length; k++) {
		    s[i * sInit.length + k] = initModel.getValues(sInit[k]);
		    y[i * sInit.length + k] = initModel.getValues(yInit[k]);
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
	return new Solution(s, x, y, store);
    }

    private static Solution reOptimize(Solution init, Store store, HAPSA.Objective obj, double param) {
	System.out.println("Initiating " + obj + " re-optimization procedure...");

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
	Solution.Objective solObj;
	switch (obj) {
	case AVA:
	    solObj = Solution.Objective.AVA;
	    break;
	case HLUR:
	    solObj = Solution.Objective.HLUR;
	    break;
	default:
	    solObj = Solution.Objective.VIS;
	    break;
	}
	double objective = incumbent.getObjective(solObj, param);
	System.out.println("Objective = " + objective);

	class ObjectiveComparator implements Comparator<Shelf> {
	    @Override
	    public int compare(Shelf shelf1, Shelf shelf2) {
		double obj1 = incumbent.getShelfObjective(shelf1, solObj, param);
		double obj2 = incumbent.getShelfObjective(shelf2, solObj, param);
		return Double.compare(obj1, obj2);
	    }
	}

	int loops = 0;
	long startTime = System.currentTimeMillis() / 1000;
	long time = System.currentTimeMillis() / 1000;
	while (((upperBound - objective) / upperBound <= STOP_GAP) && (loops < STOP_LOOPS)
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
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    int[] shelfIndices = new int[N_REOPT];
		    int nk = partialStore.getSegments().size();
		    int nj = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][nj];
		    double[][] partialS = new double[nk][nj];
		    double[][] partialY = new double[nk][nj];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nk; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double newObjective = incumbent.getObjective(solObj, param);
		    if (newObjective < objective) {
			loops++;
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
	    time = System.currentTimeMillis();
	}
	return incumbent;
    }

    private static Solution reOptimizeAPSA(Solution init, Store store) {
	System.out.println("Initiating APSA re-optimization procedure...");

	// Solve continuous relaxation of the model.
	double upperBound = Double.MAX_VALUE;
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setObjective(HAPSA.Objective.APSA);
	    continuous.relax();
	    boolean feasible = continuous.solve();
	    if (!feasible) {
		throw new IllegalStateException("No feasible upper bound found.");
	    }
	    upperBound = continuous.getObjValue();
	    // TODO: checker
	    double[][] cS = new double[store.getSegments().size()][store.getProducts().size()];
	    double[][] cX = new double[store.getShelves().size()][store.getProducts().size()];
	    double[][] cY = new double[store.getSegments().size()][store.getProducts().size()];
	    double[] sumS = new double[store.getSegments().size()];
	    double[] sumX = new double[store.getShelves().size()];
	    double[] sumY = new double[store.getSegments().size()];

	    for (int i = 0; i < store.getShelves().size(); i++) {
		cX[i] = continuous.getValues(continuous.getX()[i]);
		for (int j = 0; j < store.getProducts().size(); j++) {
		    sumX[i] += cX[i][j];
		}
	    }
	    for (int k = 0; k < store.getSegments().size(); k++) {
		cS[k] = continuous.getValues(continuous.getS()[k]);
		cY[k] = continuous.getValues(continuous.getY()[k]);
		for (int j = 0; j < store.getProducts().size(); j++) {
		    sumS[k] += cS[k][j];
		    sumY[k] += cY[k][j];
		}
	    }
	    System.out.println("Upper bound: " + upperBound + " (" + continuous.getStatus() + ")");
	} catch (IloException e) {
	    System.out.println("Continuous relaxation could not be solved in Main.reOptimize(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}

	Solution incumbent = init;
	double objective = incumbent.getObjectiveAPSA();
	System.out.println("Objective = " + objective);

	class ObjectiveComparator implements Comparator<Shelf> {
	    @Override
	    public int compare(Shelf shelf1, Shelf shelf2) {
		double obj1 = incumbent.getShelfObjectiveAPSA(shelf1);
		double obj2 = incumbent.getShelfObjectiveAPSA(shelf2);
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

	    // TODO: checker
	    ArrayList<Double> objectives = new ArrayList<Double>(shelves.size());
	    for (Shelf shelf : shelves) {
		objectives.add(incumbent.getShelfObjectiveAPSA(shelf));
	    }
	    double[] incS = new double[incumbent.getS().length];
	    double[] incX = new double[incumbent.getX().length];
	    double[] incY = new double[incumbent.getY().length];

	    for (int j = 0; j < incumbent.getS()[0].length; j++) {
		for (int i = 0; i < incX.length; i++) {
		    incX[i] += incumbent.getX()[i][j];
		}
		for (int k = 0; k < incS.length; k++) {
		    incS[k] += incumbent.getS()[k][j];
		    incY[k] += incumbent.getY()[k][j];
		}
	    }

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

		// TODO: checker
		ArrayList<Double> selectedObjectives = new ArrayList<Double>(selected.size());
		for (Shelf shelf : selected) {
		    selectedObjectives.add(incumbent.getShelfObjectiveAPSA(shelf));
		}

		ArrayList<Integer> consideredProducts = incumbent.getConsidered(selected);
		Store partialStore = incumbent.partialStore(selected, consideredProducts);
		try (HAPSA model = new HAPSA(partialStore)) {
		    model.initializePartial(incumbent, consideredProducts);
		    model.setObjective(HAPSA.Objective.APSA);
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    int[] shelfIndices = new int[N_REOPT];
		    int nk = partialStore.getSegments().size();
		    int nj = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][nj];
		    double[][] partialS = new double[nk][nj];
		    double[][] partialY = new double[nk][nj];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nk; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    double partialObj = model.getObjValue();
		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double partialObj2 = 0;
		    for (Shelf shelf : selected) {
			double shelfobj = incumbent.getShelfObjectiveAPSA(shelf);
			partialObj2 += incumbent.getShelfObjectiveAPSA(shelf);
		    }
		    double newObjective = incumbent.getObjectiveAPSA();
		    if (newObjective < objective) {
			loops++;
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
	    time = System.currentTimeMillis();
	}
	return incumbent;
    }

    private static Solution reOptimizeHAPSA(Solution init, Store store, double gamma, double theta) {
	System.out.println("Initiating HAPSA re-optimization procedure...");

	// Solve continuous relaxation of the model.
	double upperBound = Double.MAX_VALUE;
	try (HAPSA continuous = new HAPSA(store)) {
	    continuous.setGamma(gamma);
	    continuous.setTheta(theta);
	    continuous.setObjective(HAPSA.Objective.HAPSA);
	    continuous.relax();
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
	double objective = incumbent.getObjectiveHAPSA(gamma, theta);
	System.out.println("Objective = " + objective);

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
	while (((upperBound - objective) / upperBound <= STOP_GAP) && (loops < STOP_LOOPS)
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
		    model.setObjective(HAPSA.Objective.HAPSA);
		    boolean feasible = model.solve();

		    if (!feasible) {
			throw new IllegalStateException("No feasible solution found.");
		    }

		    System.out.println("Feasible solution found! Type: " + model.getStatus() + "\n");

		    int[] shelfIndices = new int[N_REOPT];
		    int nk = partialStore.getSegments().size();
		    int nj = partialStore.getProducts().size();
		    double[][] partialX = new double[N_REOPT][nj];
		    double[][] partialS = new double[nk][nj];
		    double[][] partialY = new double[nk][nj];

		    for (int i = 0; i < N_REOPT; i++) {
			shelfIndices[i] = store.getShelves().indexOf(selected.get(i));
			partialX[i] = model.getValues(model.getX()[i]);
		    }

		    for (int k = 0; k < nk; k++) {
			partialS[k] = model.getValues(model.getS()[k]);
			partialY[k] = model.getValues(model.getY()[k]);
		    }

		    incumbent.updatePartial(shelfIndices, consideredProducts, partialS, partialX, partialY);
		    double newObjective = incumbent.getObjectiveHAPSA(gamma, theta);
		    if (newObjective < objective) {
			loops++;
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
	    time = System.currentTimeMillis();
	}
	return incumbent;
    }
}
