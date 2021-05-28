/*
 * SSP.java
 * 
 * v1.0
 *
 * 13 May 2021
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

import ilog.concert.IloException;
import ilog.concert.IloIntExpr;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
 * The SSP class provides a framework for the Single Shelf Problem (SSP) model
 * used in the initialization procedure of the the optimization-based heuristic
 * approach to assortment planning and shelf space optimization. It extends the
 * {@link IloCplex} class.
 * 
 * The decision variables and constraints closely follow the notation and
 * numbering used in Flamand (2018), see:
 * 
 * T. Flamand, A. Ghoniem, M. Haouari, B. Maddah, Integrated assortment planning
 * and store-wide shelf space allocation: An optimization-based approach, Omega
 * 81 (2018) 134–149. doi:10.1016/j.omega.2017.10.006.
 *
 * @author Stefan van Berkum
 *
 */
public class SSP extends IloCplex {

    /**
     * The possible objective types: Assortment Planning and Shelf-space Allocation
     * (APSA), APSA with an availability penalty (AVA), Health-adjusted APSA
     * (HAPSA), APSA with a Healthy-Left, Unhealthy Right approach, and APSA with a
     * visibility penalty (VIS).
     *
     * @author Stefan van Berkum
     *
     */
    enum Objective {
	APSA, AVA, HAPSA, HLUR, VIS
    }

    /**
     * Automatically generated serial ID.
     */
    private static final long serialVersionUID = 4942197096268640510L;

    /**
     * The objective-specific parameter for the objectives with a visibility penalty
     * (VIS and HAPSA).
     */
    private Double gamma;

    /**
     * The objective-specific parameter for the objectives with an availability
     * penalty (AVA).
     */
    private Double lambda;

    /**
     * Decision variable q_kj, equal to one iff product category j is assigned to
     * both segments k and k+1.
     */
    private final IloIntVar[][] q;

    /**
     * Decision variable s_kj, the amount of space allocated to product category j
     * on shelf segment k.
     */
    private final IloNumVar[][] s;

    /** The shelf that is considered in this SSP model. */
    private final Shelf shelf;

    /** The store that is considered in this SSP model. */
    private final Store store;

    /**
     * The objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    private Double theta;

    /**
     * Decision variable w_j, equal to one iff product category j is assigned to
     * this shelf.
     */
    private final IloIntVar[] w;

    /**
     * Decision variable y_kj, equal to one iff product category j is assigned to
     * shelf segment k.
     */
    private final IloIntVar[][] y;

    /**
     * Constructs an SSP model.
     * 
     * @param store the single shelf store that is considered in this SSP model
     * @throws IloException if the instance could not be created
     */
    public SSP(Store store) throws IloException {
	super();
	this.store = store;
	this.shelf = store.getShelves().get(0);
	this.q = new IloIntVar[shelf.getSegments().size()][store.getProducts().size()];
	this.initBool(this.q);
	this.s = new IloNumVar[shelf.getSegments().size()][store.getProducts().size()];
	this.initNum(this.s, 0, Double.MAX_VALUE);
	this.w = new IloIntVar[store.getProducts().size()];
	this.initBool(this.w);
	this.y = new IloIntVar[shelf.getSegments().size()][store.getProducts().size()];
	this.initBool(this.y);
	this.addConstraints();
    }

    /**
     * Gets the decision variable s_kj.
     *
     * @return the decision variable s_kj
     */
    public IloNumVar[][] getS() {
	return this.s;
    }

    /**
     * Gets the decision variable w_j.
     *
     * @return the decision variable w_j
     */
    public IloIntVar[] getW() {
	return this.w;
    }

    /**
     * Gets the decision variable y_kj.
     *
     * @return the decision variable y_kj
     */
    public IloIntVar[][] getY() {
	return this.y;
    }

    /**
     * Sets the objective-specific parameter for the objectives with a visibility
     * penalty (VIS and HAPSA).
     */
    public void setGamma(double val) {
	this.gamma = val;
    }

    /**
     * Sets the objective-specific parameter for the objectives with an availability
     * penalty (AVA).
     */
    public void setLambda(double val) {
	this.lambda = val;
    }

    /**
     * Set the objective function to the user-defined type.
     * 
     * @param obj the objective function that needs to be maximized in this SSP
     *            model, one of: APSA, AVA, HAPSA, HLUR, VIS
     */
    public void setObjective(Objective obj) {
	try {
	    IloNumExpr objective = this.generateObjective(obj);
	    this.addMaximize(objective);
	} catch (IloException e) {
	    System.err.println("Objective could not be initialized.");
	    e.printStackTrace();
	    System.exit(-1);
	} catch (IllegalStateException e) {
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Sets the objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    public void setTheta(double val) {
	this.theta = val;
    }

    /**
     * Equation (4b) in Flamand (2018). This constraint ensures that the allocated
     * space for any given segment is smaller or equal to its capacity.
     * 
     * @param segment any segment
     * @param k       the index of the segment
     */
    private void add4b(Segment segment, int k) {
	try {
	    double c = segment.getCapacity();
	    this.addLe(this.sum(this.s[k]), c, "4b" + (k + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4b could not be added in SSP.add4b(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4c) in Flamand (2018). This constraint ensures that any product
     * selected on this shelf is allocated an amount of space between its minimum
     * and maximum space requirement.
     * 
     * @param sT      the transpose of the matrix of decision variables s_kj
     * @param product any product
     * @param j       the index of the product
     */
    private void add4c(IloNumVar[][] sT, Product product, int j) {
	try {
	    IloNumExpr lw = this.prod(product.getMinSpace(), this.w[j]);
	    IloNumExpr uw = this.prod(product.getMaxSpace(), this.w[j]);
	    this.addGe(this.sum(sT[j]), lw, "4ca" + (j + 1));
	    this.addLe(this.sum(sT[j]), uw, "4cb" + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4c could not be added in SSP.add4c(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4d) in Flamand (2018). This constraint ensures that the allocated
     * space for any given product on any given shelf segment is between the minimum
     * allocated space for this product, and the minimum of the capacity of this
     * segment and the maximum space requirement of this product.
     * 
     * @param segment any segment
     * @param k       the index of the segment
     * @param product any product
     * @param j       the index of the product
     */
    private void add4d(Segment segment, int k, Product product, int j) {
	try {
	    IloNumExpr psiy = this.prod(product.getMinAllocated(), this.y[k][j]);
	    double mincu = Math.min(segment.getCapacity(), product.getMaxSpace());
	    IloNumExpr mincuy = this.prod(mincu, this.y[k][j]);
	    this.addGe(this.s[k][j], psiy, "4da" + (k + 1) + (j + 1));
	    this.addLe(this.s[k][j], mincuy, "4db" + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4d could not be added in SSP.add4d(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4e) in Flamand (2018). This constraint ensures that any product
     * that is allocated to a pair of segments k1 and k3 (k1 < k2 < k3), is also
     * allocated to any segment in between (k2).
     * 
     * @param k1      the index of any segment k1
     * @param product any product
     * @param j       the index of the product
     */
    private void add4e(int k1, Product product, int j) {
	try {
	    for (int k3 = k1 + 2; k3 < this.shelf.getSegments().size(); k3++) {
		// Check whether (k1, k3, j) in R.
		double sumc = 0;
		for (int k2 = k1 + 1; k2 < k3; k2++) {
		    sumc += this.shelf.getSegments().get(k2).getCapacity();
		}
		if (sumc <= product.getMaxSpace() - 2 * product.getMinAllocated()) {
		    // (k1, k3, j) is not in R.
		    for (int k2 = k1 + 1; k2 < k3; k2++) {
			double c2 = this.shelf.getSegments().get(k2).getCapacity();
			IloIntExpr yy1 = this.sum(this.sum(this.y[k1][j], this.y[k3][j]), -1);
			IloNumExpr cyy1 = this.prod(c2, yy1);
			this.addGe(s[k2][j], cyy1, "4e" + (k1 + 1) + (k2 + 1) + (k3 + 1) + (j + 1));
		    }
		}
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 4e could not be added in SSP.add4e(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4f) in Flamand (2018). This constraint ensures that a product is
     * only assigned to a shelf segment if it is assigned to this shelf.
     * 
     * @param k the index of any segment
     * @param j the index of any product
     */
    private void add4f(int k, int j) {
	try {
	    this.addLe(this.y[k][j], this.w[j], "4f" + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4f could not be added in SSP.add4f(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4g) in Flamand (2018). This constraint ensures that a product is
     * only assigned to this shelf if it is assigned to at least one shelf segment.
     * 
     * @param yT the transpose of the matrix of decision variables y_kj
     * @param j  the index of any product
     */
    private void add4g(IloIntVar[][] yT, int j) {
	try {
	    this.addLe(this.w[j], this.sum(yT[j]), "4g" + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4g could not be added in SSP.add4g(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4h) in Flamand (2018). This constraint ensures that q_kj is equal
     * to one iff product j is assigned to both segments k and k + 1, and zero
     * otherwise.
     * 
     * @param k the index of any segment
     * @param j the index of any product
     */
    private void add4h(int k, int j) {
	try {
	    IloIntExpr yy1 = this.sum(this.sum(this.y[k][j], this.y[k + 1][j]), -1);
	    this.addGe(this.q[k][j], yy1, "4h" + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4h could not be added in SSP.add4h(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4i) in Flamand (2018). This constraint ensures that only one
     * product category runs over any two consecutive shelf segments.
     * 
     * @param k the index of any segment
     */
    private void add4i(int k) {
	try {
	    this.addLe(this.sum(this.q[k]), 1, "4i" + (k + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 4i could not be added in SSP.add4i(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4j) in Flamand (2018). This constraint ensures that any two
     * products that have allocation disaffinity, are not placed on the same shelf.
     */
    private void add4j() {
	try {
	    for (SymmetricPair pair : this.store.getAllocationDisaffinity()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		this.addLe(this.sum(this.w[j], this.w[jprime]), 1, "4j" + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 4j could not be added in SSP.add4j(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4k) in Flamand (2018). This constraint ensures that any two
     * products that have symmetric assortment affinity, are either selected
     * together or neither, and placed on the same shelf.
     */
    private void add4k() {
	try {
	    for (SymmetricPair pair : this.store.getSymmetricAssortment()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		IloIntExpr ww = this.diff(this.w[j], this.w[jprime]);
		this.addEq(ww, 0, "4k" + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 4k could not be added in SSP.add4k(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4l) in Flamand (2018). This constraint ensures that for any two
     * products that have asymmetric assortment affinity: if the first product is
     * selected, the second product is selected as well, and placed on the same
     * shelf.
     */
    private void add4l() {
	try {
	    for (AsymmetricPair pair : this.store.getAsymmetricAssortment()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		this.addLe(this.w[j], this.w[jprime], "4l" + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 4l could not be added in SSP.add4l(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (4m) in Flamand (2018). This constraint ensures that any product
     * that cannot cover all intermediate sections between segments k1 and k3 (k1 <
     * k2 < k3), is not allocated to both segments k1 and k3.
     */
    private void add4m(int k1, Product product, int j) {
	try {
	    for (int k3 = k1 + 2; k3 < this.shelf.getSegments().size(); k3++) {
		// Check whether (k1, k3, j) in R.
		double sumc = 0;
		for (int k2 = k1 + 1; k2 < k3; k2++) {
		    sumc += this.shelf.getSegments().get(k2).getCapacity();
		}
		if (sumc > product.getMaxSpace() - 2 * product.getMinAllocated()) {
		    // (k1, k3, j) is in R.
		    IloIntExpr yy = this.sum(this.y[k1][j], this.y[k3][j]);
		    this.addLe(yy, 1, "4m" + (k1 + 1) + (k3 + 1) + (j + 1));
		}
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 4j could not be added in SSP.add4j(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Adds constraints for this SSP model.
     */
    private void addConstraints() {
	// Get transpose of all matrix variables.
	IloNumVar[][] sT = Utils.getTranspose(this.s);
	IloIntVar[][] yT = Utils.getTranspose(this.y);

	// All constraints with: (for each segment).
	for (int k = 0; k < this.shelf.getSegments().size(); k++) {
	    Segment segment = this.shelf.getSegments().get(k);

	    this.add4b(segment, k);

	    if (k < this.shelf.getSegments().size() - 1) {
		this.add4i(k);
	    }

	    // All constraints with: (for each segment, for each product).
	    for (int j = 0; j < this.store.getProducts().size(); j++) {
		Product product = this.store.getProducts().get(j);

		this.add4d(segment, k, product, j);
		this.add4e(k, product, j);
		this.add4f(k, j);
		this.add4m(k, product, j);

		if (k < this.shelf.getSegments().size() - 1) {
		    this.add4h(k, j);
		}
	    }
	}

	// All constraints with: (for each product).
	for (int j = 0; j < this.store.getProducts().size(); j++) {
	    Product product = this.store.getProducts().get(j);

	    this.add4c(sT, product, j);
	    this.add4g(yT, j);
	}

	// Affinity relationship constraints.
	this.add4j();
	this.add4k();
	this.add4l();
    }

    /**
     * Generates the objective value based on a user-defined type.
     * 
     * @param obj the objective value type
     * @return the objective value to be maximized
     * @throws IllegalStateException if the required parameters are not initialized
     */
    private IloNumExpr generateObjective(Objective obj) throws IllegalStateException {
	try {
	    IloNumExpr minuend = this.constant(0); // The base objective (profit).
	    IloNumExpr subtrahend = this.constant(0); // The objective-specific subtrahend (zero for APSA).
	    for (int k = 0; k < this.shelf.getSegments().size(); k++) {
		Segment segment = this.shelf.getSegments().get(k);

		for (int j = 0; j < this.store.getProducts().size(); j++) {
		    Product product = this.store.getProducts().get(j);

		    double fc = segment.getAttractiveness() / segment.getCapacity();
		    IloNumExpr fsc = this.prod(fc, this.s[k][j]);
		    IloNumExpr phifsc = this.prod(product.getMaxProfit(), fsc);
		    minuend = this.sum(minuend, phifsc);

		    switch (obj) {
		    case AVA:
			if (this.lambda == null) {
			    throw new IllegalStateException("Lambda has not been initialized.");
			}

			double lambdah = this.lambda / product.getHealthScore();
			subtrahend = this.sum(subtrahend, this.prod(lambdah, this.y[k][j]));
			break;
		    case HAPSA:
			if (this.gamma == null || this.theta == null) {
			    throw new IllegalStateException("Gamma or theta has not been initialized.");
			}

			// Visibility penalty.
			double gammah = this.gamma / product.getHealthScore();
			subtrahend = this.sum(subtrahend, this.prod(gammah, fsc));

			// Healthy-left, unhealthy-right approach.
			double thetah = this.theta * product.getHealthScore();
			minuend = this.sum(minuend, this.prod(thetah, this.y[k][j]));
			break;
		    case HLUR:
			if (this.theta == null) {
			    throw new IllegalStateException("Theta has not been initialized.");
			}

			double thetah2 = this.theta * product.getHealthScore();
			minuend = this.sum(minuend, this.prod(thetah2, this.y[k][j]));
			break;
		    case VIS:
			if (this.gamma == null) {
			    throw new IllegalStateException("Gamma has not been initialized.");
			}

			double gammah2 = this.gamma / product.getHealthScore();
			subtrahend = this.sum(subtrahend, this.prod(gammah2, fsc));
			break;
		    default:
			break;
		    }
		}
		if (obj == Objective.HAPSA || obj == Objective.HLUR) {
		    IloIntExpr sumsumy2h2 = this.constant(0);
		    int nh = this.shelf.getHorizontal();

		    if ((k + 1) % nh != 0) {
			// Shelf segment k is not the rightmost segment.

			for (int k2 = k + 1; k2 <= k + (nh - (k + 1) % nh); k2++) {
			    // Shelf segments to the right of shelf segment k.

			    for (int j2 = 0; j2 < this.store.getProducts().size(); j2++) {
				Product product2 = this.store.getProducts().get(j2);
				IloIntExpr y2h2 = this.prod(this.y[k2][j2], product2.getHealthScore());
				sumsumy2h2 = this.sum(sumsumy2h2, y2h2);
			    }
			}
			IloNumExpr thetasumsumy2h2 = this.prod(this.theta, sumsumy2h2);
			subtrahend = this.sum(subtrahend, thetasumsumy2h2);
		    }
		}
	    }
	    return this.diff(minuend, subtrahend);
	} catch (IloException e) {
	    System.err.println("Objective could not be created in SSP.getObjective.");
	    e.printStackTrace();
	    System.exit(-1);
	    return null;
	}
    }

    /**
     * Initializes a vector of boolean decision variables.
     * 
     * @param var the vector of boolean variables
     */
    private void initBool(IloIntVar[] var) {
	for (int i = 0; i < var.length; i++) {
	    try {
		var[i] = this.boolVar();
	    } catch (IloException e) {
		System.err.println("A boolean variable could not be created in SSP.initBool(...).");
		e.printStackTrace();
		System.exit(-1);
	    }
	}
    }

    /**
     * Initializes a matrix of boolean decision variables.
     * 
     * @param var the matrix of boolean variables
     */
    private void initBool(IloIntVar[][] var) {
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		try {
		    var[i][j] = this.boolVar();
		} catch (IloException e) {
		    System.err.println("A boolean variable could not be created in SSP.initBool(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	}
    }

    /**
     * Initializes a matrix of number decision variables.
     * 
     * @param var a matrix of number variables
     * @param lb  the lower bound for each variable
     * @param ub  the upper bound for each variable
     */
    private void initNum(IloNumVar[][] var, double lb, double ub) {
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		try {
		    var[i][j] = this.numVar(lb, ub);
		} catch (IloException e) {
		    System.err.println("A number variable could not be created in SSP.initNum(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	}
    }
}
