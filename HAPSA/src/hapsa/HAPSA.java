/*
 * HAPSA.java
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

import java.util.ArrayList;

import ilog.concert.IloConversion;
import ilog.concert.IloException;
import ilog.concert.IloIntExpr;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;

/**
 * The HAPSA class provides a framework for the Health-adjusted Assortment
 * Planning and Shelf-space Allocation (HAPSA) model used in the MIP-based
 * re-optimization procedure of the optimization-based heuristic approach to
 * assortment planning and shelf space optimization. It extends the
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
public class HAPSA extends IloCplex {

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
    private static final long serialVersionUID = -2862092852120058808L;

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
    private IloIntVar[][] q;

    /**
     * Decision variable s_kj, the amount of space allocated to product category j
     * on shelf segment k.
     */
    private IloNumVar[][] s;

    /** The store that is considered in this HAPSA model. */
    private Store store;

    /**
     * The objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    private Double theta;

    /**
     * Decision variable x_ij, equal to one iff product category j is assigned to
     * shelf i.
     */
    private IloIntVar[][] x;

    /**
     * Decision variable y_kj, equal to one iff product category j is assigned to
     * shelf segment k.
     */
    private IloIntVar[][] y;

    /**
     * Decision variable z_jj', equal to one iff product categories j and j' are
     * selected in the assortment simultaneously.
     */
    private IloIntVar[][] z;

    /**
     * Constructs a HAPSA model.
     * 
     * @param store the store that is considered in this HAPSA model
     * @throws IloException if the instance could not be created
     */
    public HAPSA(Store store) throws IloException {
	super();
	this.store = store;
	this.q = new IloIntVar[store.getSegments().size()][store.getProducts().size()];
	initBool(this.q);
	this.s = new IloNumVar[store.getSegments().size()][store.getProducts().size()];
	initNum(this.s, 0, Double.MAX_VALUE);
	this.x = new IloIntVar[store.getShelves().size()][store.getProducts().size()];
	initBool(this.x);
	this.y = new IloIntVar[store.getSegments().size()][store.getProducts().size()];
	initBool(this.y);
	this.z = new IloIntVar[store.getProducts().size()][store.getProducts().size()];
	initBool(this.z);
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
     * Gets the decision variable x_ij.
     *
     * @return the decision variable x_ij
     */
    public IloIntVar[][] getX() {
	return this.x;
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
     * Initializes the HAPSA model decision variables s_kj, x_ij, and y_kj for a
     * partial model.
     * 
     * @param incumbent the incumbent solution
     * @param products  the list of considered product indices
     */
    public void initializePartial(Solution incumbent, ArrayList<Integer> products) {
	try {
	    Store store = incumbent.getStore();
	    double[][] s = incumbent.getS();
	    double[][] x = incumbent.getX();
	    double[][] y = incumbent.getY();
	    int ni = this.store.getShelves().size();
	    int nk = this.store.getSegments().size();
	    int nj = this.store.getProducts().size();
	    double[] initParam = new double[nj * (2 * nk + ni)]; // Concatenation of s, x, and y.

	    // Collect partial decision variables.
	    int index = 0;
	    for (int seg = 0; seg < nk; seg++) {
		Segment segment = this.store.getSegments().get(seg);
		int k = store.getSegments().indexOf(segment);
		for (int pr = 0; pr < nj; pr++) {
		    int j = products.get(pr);
		    initParam[index++] = s[k][j];
		}
	    }
	    for (int sh = 0; sh < ni; sh++) {
		Shelf shelf = this.store.getShelves().get(sh);
		int i = store.getShelves().indexOf(shelf);
		for (int pr = 0; pr < nj; pr++) {
		    int j = products.get(pr);
		    initParam[index++] = x[i][j];
		}
	    }
	    for (int seg = 0; seg < nk; seg++) {
		Segment segment = this.store.getSegments().get(seg);
		int k = store.getSegments().indexOf(segment);
		for (int pr = 0; pr < nj; pr++) {
		    int j = products.get(pr);
		    initParam[index++] = y[k][j];
		}
	    }

	    // Initialize partial model.
	    IloNumVar[] modelParam = Utils.concatenate(this.s, this.x, this.y);
	    this.addMIPStart(modelParam, initParam);
	} catch (

	IloException e) {
	    System.err.println("Parameters could not be initialized in HAPSA.initialize(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Instructs the model to chill out. Relaxes all integer variables to be
     * continuous.
     */
    public void relax() {
	try {
	    for (int i = 0; i < this.store.getShelves().size(); i++) {
		IloConversion xConv = this.conversion(x[i], IloNumVarType.Float);
		this.add(xConv);
	    }
	    for (int k = 0; k < this.store.getSegments().size(); k++) {
		IloConversion qConv = this.conversion(q[k], IloNumVarType.Float);
		this.add(qConv);
		IloConversion yConv = this.conversion(y[k], IloNumVarType.Float);
		this.add(yConv);
	    }
	    for (int j = 0; j < this.store.getProducts().size(); j++) {
		IloConversion zConv = this.conversion(z[j], IloNumVarType.Float);
		this.add(zConv);
	    }
	} catch (IloException e) {
	    System.err.println("Model could not relax in HAPSA.relax().");
	    e.printStackTrace();
	    System.exit(-1);
	}
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
     * Sets the decision variable s_kj.
     *
     * @param s the decision variable s_kj to set
     */
    public void setS(IloNumVar[][] s) {
	this.s = s;
    }

    /**
     * Sets the objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    public void setTheta(double val) {
	this.theta = val;
    }

    /**
     * Sets the decision variable x_ij.
     *
     * @param x the decision variable x_ij to set
     */
    public void setX(IloIntVar[][] x) {
	this.x = x;
    }

    /**
     * Sets the decision variable y_kj.
     *
     * @param y the decision variable y_kj to set
     */
    public void setY(IloIntVar[][] y) {
	this.y = y;
    }

    /**
     * Equation (1b) in Flamand (2018). This constraint ensures that each product is
     * assigned to at most one shelf.
     * 
     * @param xT the transpose of the matrix of decision variables x_ij
     * @param j  the index of any product
     */
    private void add1b(IloIntVar[][] xT, int j) {
	try {
	    this.addLe(this.sum(xT[j]), 1, "1b" + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1b could not be added in HAPSA.add1b(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1c) in Flamand (2018). This constraint ensures that the allocated
     * space on any given shelf does not exceed its capacity.
     * 
     * @param segment any segment
     * @param k       the index of the segment
     */
    private void add1c(Segment segment, int k) {
	try {
	    this.addLe(this.sum(this.s[k]), segment.getCapacity(), "1c" + (k + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1c could not be added in HAPSA.add1c(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1d) in Flamand (2018). This constraint ensures that the allocated
     * space for each product is between its minimum and maximum space requirement.
     * 
     * @param sT      the transpose of the matrix of decision variables s_kj
     * @param xT      the transpose of the matrix of decision variables x_ij
     * @param product any product
     * @param j       the index of the product
     */
    private void add1d(IloNumVar[][] sT, IloIntVar[][] xT, Product product, int j) {
	try {
	    IloNumExpr lsumx = this.prod(product.getMinSpace(), this.sum(xT[j]));
	    IloNumExpr usumx = this.prod(product.getMaxSpace(), this.sum(xT[j]));
	    this.addGe(this.sum(sT[j]), lsumx, "1da" + (j + 1));
	    this.addLe(this.sum(sT[j]), usumx, "1db" + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1d could not be added in HAPSA.add1d(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1e) in Flamand (2018). This constraint ensures that the allocated
     * space for any given product on any given shelf segment is between the minimum
     * allocated space for this product, and the minimum of the capacity of this
     * segment and the maximum space requirement of this product.
     * 
     * @param segment any segment
     * @param k       the index of the segment
     * @param product any product
     * @param j       the index of the product
     */
    private void add1e(Segment segment, int k, Product product, int j) {
	try {
	    IloNumExpr psiy = this.prod(product.getMinAllocated(), this.y[k][j]);
	    double mincu = Math.min(segment.getCapacity(), product.getMaxSpace());
	    IloNumExpr mincuy = this.prod(mincu, this.y[k][j]);
	    this.addGe(this.s[k][j], psiy, "1ea" + (k + 1) + (j + 1));
	    this.addLe(this.s[k][j], mincuy, "1eb" + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1e could not be added in HAPSA.add1e(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1g) in Flamand (2018). This constraint ensures that any product is
     * only allocated to a shelf segment if it is assigned to the corresponding
     * shelf.
     * 
     * @param shelf any shelf
     * @param i     the index of any shelf
     * @param k     the index of any segment on the shelf
     * @param j     the index of any product
     */
    private void add1g(Shelf shelf, int i, int k, int j) {
	try {
	    int startK = i * shelf.getSegments().size();
	    this.addLe(this.y[startK + k][j], this.x[i][j], "1g" + (i + 1) + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1g could not be added in HAPSA.add1g(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1h) in Flamand (2018). This constraint ensures that any product is
     * placed on at least one shelf segment if it is allocated to the corresponding
     * shelf.
     * 
     * @param yT the transpose of the matrix of decision variables y_kj
     * @param i  the index of any shelf
     * @param j  the index of any product
     */
    private void add1h(IloIntVar[][] yT, Shelf shelf, int i, int j) {
	try {
	    int startIndex = i * shelf.getSegments().size();
	    int nSegments = shelf.getSegments().size();
	    this.addLe(this.x[i][j], this.sum(yT[j], startIndex, nSegments), "1h" + (i + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1h could not be added in HAPSA.add1h(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1i) in Flamand (2018). This constraint ensures that for any shelf
     * q_kj is equal to one iff product j is assigned to both segments k and k + 1,
     * and zero otherwise.
     * 
     * @param shelf any shelf
     * @param i     the index of any shelf
     * @param k     the index of any segment on the shelf
     * @param j     the index of any product
     */
    private void add1i(Shelf shelf, int i, int k, int j) {
	try {
	    int startK = i * shelf.getSegments().size();
	    IloIntExpr yy1 = this.sum(this.sum(this.y[startK + k][j], this.y[startK + (k + 1)][j]), -1);
	    this.addGe(this.q[startK + k][j], yy1, "1i" + (i + 1) + (k + 1) + (j + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1i could not be added in HAPSA.add1i(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1j) in Flamand (2018). This constraint ensures that only one
     * product category runs over any two consecutive shelf segments.
     * 
     * @param shelf any shelf
     * @param i     the index of any shelf
     * @param k     the index of any segment on the shelf
     */
    private void add1j(Shelf shelf, int i, int k) {
	try {
	    int startK = i * shelf.getSegments().size();
	    this.addLe(this.sum(this.q[startK + k]), 1, "1j" + (i + 1) + (k + 1));
	} catch (IloException e) {
	    System.err.println("Constraint 1j could not be added in HAPSA.add1j(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1k) in Flamand (2018). This constraint ensures that any two
     * products that have allocation disaffinity, are not placed on the same shelf.
     * 
     * @param i the index of any shelf
     */
    private void add1k(int i) {
	try {
	    for (SymmetricPair pair : this.store.getAllocationDisaffinity()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		this.addLe(this.sum(this.x[i][j], this.x[i][jprime]), 1, "1k" + (i + 1) + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 1k could not be added in HAPSA.add1k(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1l) in Flamand (2018). This constraint ensures that any two
     * products that have symmetric assortment affinity, are either selected
     * together or neither, and placed on the same shelf.
     * 
     * @param i the index of any shelf
     */
    private void add1l(int i) {
	try {
	    for (SymmetricPair pair : this.store.getSymmetricAssortment()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		IloIntExpr xx = this.diff(this.x[i][j], this.x[i][jprime]);
		this.addEq(xx, 0, "1l" + (i + 1) + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 1l could not be added in HAPSA.add1l(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (1m) in Flamand (2018). This constraint ensures that for any two
     * products that have asymmetric assortment affinity: if the first product is
     * selected, the second product is selected as well, and placed on the same
     * shelf.
     * 
     * @param i the index of any shelf
     */
    private void add1m(int i) {
	try {
	    for (AsymmetricPair pair : this.store.getAsymmetricAssortment()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();
		this.addLe(this.x[i][j], this.x[i][jprime], "1m" + (i + 1) + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 1m could not be added in HAPSA.add1m(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equations (1n) and (1o) in Flamand (2018). These constraints ensure that any
     * two products that have allocation affinity are placed on the same shelf.
     * 
     * @param i the index of any shelf
     */
    private void add1no(int i) {
	try {
	    for (SymmetricPair pair : this.store.getAllocationAffinity()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();

		// Equation (1n).
		IloIntExpr xx = this.diff(this.x[i][j], this.x[i][jprime]);
		this.addLe(xx, this.diff(1, this.z[j][jprime]), "1n" + (i + 1) + j + jprime);

		// Equation (1o).
		this.addGe(xx, this.sum(-1, this.z[j][jprime]), "1o" + (i + 1) + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 1n and 1o could not be added in HAPSA.add1no(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equations (1p-1r) in Flamand (2018). These constraints ensure that any two
     * products that have allocation affinity are placed on the same shelf.
     * 
     * @param xT the transpose of the matrix of decision variables x_ij
     */
    private void add1pqr(IloIntVar[][] xT) {
	try {
	    for (SymmetricPair pair : this.store.getAllocationAffinity()) {
		int j = pair.getIndex1();
		int jprime = pair.getIndex2();

		// Equation (1p).
		this.addLe(this.z[j][jprime], this.sum(xT[j]), "1p" + j);

		// Equation (1q).
		this.addLe(this.z[j][jprime], this.sum(xT[jprime]), "1q" + jprime);

		// Equation (1r).
		IloIntExpr sumxsumx = this.sum(this.sum(xT[j]), this.sum(xT[jprime]));
		IloIntExpr sumxsumx1 = this.sum(sumxsumx, -1);
		this.addGe(this.z[j][jprime], sumxsumx1, "1r" + j + jprime);
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 1p-1r could not be added in HAPSA.add1pqr(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (2) in Flamand (2018). This constraint ensures that any product that
     * cannot cover all intermediate sections between segments k1 and k3 (k1 < k2 <
     * k3) on any given shelf, is not allocated to both segments k1 and k3.
     * 
     * @param shelf   any shelf
     * @param i       the index of the shelf
     * @param k1      the index of any segment on the shelf
     * @param product any product
     * @param j       the index of the product
     */
    private void add2(Shelf shelf, int i, int k1, Product product, int j) {
	try {
	    int startK = i * shelf.getSegments().size();
	    for (int k3 = k1 + 2; k3 < shelf.getSegments().size(); k3++) {
		// Check whether (k1, k3, j) in R.
		double sumc = 0;
		for (int k2 = k1 + 1; k2 < k3; k2++) {
		    sumc += shelf.getSegments().get(k2).getCapacity();
		}
		if (sumc > product.getMaxSpace() - 2 * product.getMinAllocated()) {
		    // (k1, k3, j) is in R.
		    IloIntExpr yy = this.sum(this.y[startK + k1][j], this.y[startK + k3][j]);
		    this.addLe(yy, 1, "2" + (i + 1) + (k1 + 1) + (k3 + 1) + (j + 1));
		}
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 2 could not be added in HAPSA.add2(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Equation (3) in Flamand (2018). This constraint ensures that any product that
     * is allocated to a pair of segments k1 and k3 (k1 < k2 < k3) on any given
     * shelf, is also allocated to any segment in between (k2).
     * 
     * @param shelf   any shelf
     * @param i       the index of the shelf
     * @param k1      the index of any segment on the shelf
     * @param product any product
     * @param j       the index of the product
     */
    private void add3(Shelf shelf, int i, int k1, Product product, int j) {
	try {
	    int startK = i * shelf.getSegments().size();
	    for (int k3 = k1 + 2; k3 < shelf.getSegments().size(); k3++) {
		// Check whether (k1, k3, j) in R.
		double sumc = 0;
		for (int k2 = k1 + 1; k2 < k3; k2++) {
		    sumc += shelf.getSegments().get(k2).getCapacity();
		}
		if (sumc <= product.getMaxSpace() - 2 * product.getMinAllocated()) {
		    // (k1, k3, j) is not in R.
		    for (int k2 = k1 + 1; k2 < k3; k2++) {
			double c2 = shelf.getSegments().get(k2).getCapacity();
			IloIntExpr yy1 = this.sum(this.sum(this.y[startK + k1][j], this.y[startK + k3][j]), -1);
			IloNumExpr cyy1 = this.prod(c2, yy1);
			this.addGe(s[startK + k2][j], cyy1, "3" + (i + 1) + (k1 + 1) + (k2 + 1) + (k3 + 1) + (j + 1));
		    }
		}
	    }
	} catch (IloException e) {
	    System.err.println("Constraint 3 could not be added in HAPSA.add3(...).");
	    e.printStackTrace();
	    System.exit(-1);
	}
    }

    /**
     * Adds constraints for this HAPSA model.
     */
    private void addConstraints() {
	IloNumVar[][] sT = Utils.getTranspose(this.s);
	IloIntVar[][] xT = Utils.getTranspose(this.x);
	IloIntVar[][] yT = Utils.getTranspose(this.y);

	this.add1pqr(xT);

	// All constraints with: (for each shelf).
	for (int i = 0; i < this.store.getShelves().size(); i++) {
	    Shelf shelf = this.store.getShelves().get(i);

	    this.add1k(i);
	    this.add1l(i);
	    this.add1m(i);
	    this.add1no(i);

	    // All constraints with: (for each shelf, for each segment).
	    for (int k = 0; k < shelf.getSegments().size(); k++) {
		if (k < shelf.getSegments().size() - 1) {
		    this.add1j(shelf, i, k);
		}

		// All constraints with: (for each shelf, for each segment, for each product).
		for (int j = 0; j < this.store.getProducts().size(); j++) {
		    Product product = this.store.getProducts().get(j);

		    this.add1g(shelf, i, k, j);
		    this.add2(shelf, i, k, product, j);
		    this.add3(shelf, i, k, product, j);

		    if (k < shelf.getSegments().size() - 1) {
			this.add1i(shelf, i, k, j);
		    }
		}
	    }

	    // All constraints with: (for each shelf, for each product).
	    for (int j = 0; j < this.store.getProducts().size(); j++) {
		this.add1h(yT, shelf, i, j);
	    }
	}

	// All constraints with: (for each segment).
	for (int k = 0; k < this.store.getSegments().size(); k++) {
	    Segment segment = this.store.getSegments().get(k);

	    this.add1c(segment, k);

	    // All constraints with: (for each segment, for each product).
	    for (int j = 0; j < this.store.getProducts().size(); j++) {
		Product product = this.store.getProducts().get(j);

		this.add1e(segment, k, product, j);
	    }
	}

	// All constraints with: (for each product).
	for (int j = 0; j < this.store.getProducts().size(); j++) {
	    Product product = this.store.getProducts().get(j);

	    this.add1b(xT, j);
	    this.add1d(sT, xT, product, j);
	}
    }

    /**
     * Generates the objective value based on a user-defined type.
     * 
     * @param obj the objective value type, one of: APSA, AVA, HAPSA, HLUR, VIS
     * @return the objective value to be maximized
     * @throws IllegalStateException if the required parameters are not initialized
     */
    private IloNumExpr generateObjective(Objective obj) throws IllegalStateException {
	try {
	    IloNumExpr minuend = this.constant(0); // The base objective (profit).
	    IloNumExpr subtrahend = this.constant(0); // The objective-specific subtrahend (zero for APSA).
	    for (int i = 0; i < this.store.getShelves().size(); i++) {
		Shelf shelf = this.store.getShelves().get(i);
		int startK = i * shelf.getSegments().size();

		for (int k = 0; k < shelf.getSegments().size(); k++) {
		    Segment segment = shelf.getSegments().get(k);

		    for (int j = 0; j < this.store.getProducts().size(); j++) {
			Product product = this.store.getProducts().get(j);

			double fc = segment.getAttractiveness() / segment.getCapacity();
			IloNumExpr fsc = this.prod(fc, this.s[startK + k][j]);
			IloNumExpr phifsc = this.prod(product.getMaxProfit(), fsc);
			minuend = this.sum(minuend, phifsc);

			switch (obj) {
			case AVA:
			    if (this.lambda == null) {
				throw new IllegalStateException("Lambda has not been initialized.");
			    }

			    double lambdah = this.lambda / product.getHealthScore();
			    subtrahend = this.sum(subtrahend, this.prod(lambdah, this.y[startK + k][j]));
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
			    minuend = this.sum(minuend, this.prod(thetah, this.y[startK + k][j]));
			    break;
			case HLUR:
			    if (this.theta == null) {
				throw new IllegalStateException("Theta has not been initialized.");
			    }

			    double thetah2 = this.theta * product.getHealthScore();
			    minuend = this.sum(minuend, this.prod(thetah2, this.y[startK + k][j]));
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
			int nh = shelf.getHorizontal();

			if ((k + 1) % nh != 0) {
			    // Shelf segment k is not the rightmost segment.

			    for (int k2 = k + 1; k2 <= k + (nh - (k + 1) % nh); k2++) {
				// Shelf segments to the right of shelf segment k.

				for (int j2 = 0; j2 < this.store.getProducts().size(); j2++) {
				    Product product2 = this.store.getProducts().get(j2);
				    IloIntExpr y2h2 = this.prod(this.y[startK + k2][j2], product2.getHealthScore());
				    sumsumy2h2 = this.sum(sumsumy2h2, y2h2);
				}
			    }
			    IloNumExpr thetasumsumy2h2 = this.prod(this.theta, sumsumy2h2);
			    subtrahend = this.sum(subtrahend, thetasumsumy2h2);
			}
		    }
		}
	    }
	    return this.diff(minuend, subtrahend);
	} catch (IloException e) {
	    System.err.println("Objective could not be created in HAPSA.generateObjective.");
	    e.printStackTrace();
	    System.exit(-1);
	    return null;
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
		    System.err.println("A boolean variable could not be created in HAPSA.initBool(...).");
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
		    System.err.println("A number variable could not be created in HAPSA.initNum(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	}
    }

}
