/*
 * Model.java
 * 
 * v1.0
 *
 * 29 May 2021
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
import ilog.concert.IloIntVar;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
 * The Model class is an abstract class that defines the basic framework for a
 * model used in the optimization-based heuristic approach to assortment
 * planning and shelf space optimization. It extends the {@link IloCplex} class.
 *
 * @author Stefan van Berkum
 *
 */
public abstract class Model extends IloCplex {

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
    private static final long serialVersionUID = 1908858438240885124L;

    /**
     * The objective-specific parameter for the objectives with a visibility penalty
     * (VIS and HAPSA).
     */
    protected Double gamma;

    /**
     * The objective-specific parameter for the objectives with an availability
     * penalty (AVA).
     */
    protected Double lambda;

    /** The store that is considered in this model. */
    protected final Store store;

    /**
     * The objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    protected Double theta;

    /**
     * Constructs a basic model.
     * 
     * @param store the store considered in this model
     * @throws IloException if the instance could not be created
     */
    protected Model(Store store) throws IloException {
	this.store = store;
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
     * Sets the objective-specific parameter for the objectives with a healthy-left,
     * unhealthy-right approach (HLUR and HAPSA).
     */
    public void setTheta(double val) {
	this.theta = val;
    }

    /**
     * Initializes a vector of boolean decision variables.
     * 
     * @param var the vector of boolean variables
     */
    protected void initBool(IloIntVar[] var) {
	for (int i = 0; i < var.length; i++) {
	    try {
		var[i] = this.boolVar();
	    } catch (IloException e) {
		System.err.println("A boolean variable could not be created in Model.initBool(...).");
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
    protected void initBool(IloIntVar[][] var) {
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		try {
		    var[i][j] = this.boolVar();
		} catch (IloException e) {
		    System.err.println("A boolean variable could not be created in Model.initBool(...).");
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
    protected void initNum(IloNumVar[][] var, double lb, double ub) {
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		try {
		    var[i][j] = this.numVar(lb, ub);
		} catch (IloException e) {
		    System.err.println("A number variable could not be created in Model.initNum(...).");
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	}
    }
}
