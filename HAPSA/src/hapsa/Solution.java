/*
 * Solution.java
 * 
 * v1.0
 *
 * 19 May 2021
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
import java.util.HashSet;
import java.util.TreeSet;

/**
 * The Solution class provides a framework and methods for any store planning
 * solution.
 *
 * @author Stefan van Berkum
 *
 */
public class Solution {

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
     * Decision variable s_kj, the amount of space allocated to product category j
     * on shelf segment k.
     */
    private double[][] s;

    /* The store considered in this solution. */
    private Store store;

    /**
     * Decision variable x_ij, equal to one iff product category j is assigned to
     * shelf i.
     */
    private double[][] x;

    /**
     * Decision variable y_kj, equal to one iff product category j is assigned to
     * shelf segment k.
     */
    private double[][] y;

    /**
     * Creates a solution object, with the corresponding decision variables s_kj,
     * x_ij, and y_kj.
     * 
     * @param initS the initial decision variable s_kj
     * @param initX the initial decision variable x_ij
     * @param initY the initial decision variable y_kj
     * @param store the store
     */
    public Solution(double[][] initS, double[][] initX, double[][] initY, Store store) {
	this.s = initS;
	this.x = initX;
	this.y = initY;
	this.store = store;
    }

    /**
     * Finds and returns all products that are allocated to the given set of shelves
     * or not selected in the assortment.
     * 
     * @param shelves the set of shelves
     * @return a list of all products that are allocated to the given set of shelves
     *         or not selected in the assortment
     */
    public ArrayList<Integer> getConsidered(ArrayList<Shelf> shelves) {
	TreeSet<Integer> considered = new TreeSet<Integer>();
	// Add all products.
	for (int j = 0; j < this.store.getProducts().size(); j++) {
	    considered.add(j);
	}

	// Remove all products that are selected on other shelves.
	for (int i = 0; i < this.store.getShelves().size(); i++) {
	    Shelf shelf = this.store.getShelves().get(i);
	    if (!shelves.contains(shelf)) {
		double[] shelfX = this.x[i];

		for (int j = 0; j < this.store.getProducts().size(); j++) {
		    if (shelfX[j] == 1) {
			considered.remove(j);
		    }
		}
	    }
	}
	return new ArrayList<Integer>(considered);
    }

    /**
     * Computes the objective value of this solution for the objective functions of
     * type: AVA, HLUR, and VIS.
     * 
     * @param obj   the objective function, one of: AVA, HLUR, or VIS
     * @param param the parameter corresponding to the objective function
     * @return the objective value
     */
    public double getObjective(Objective obj, double param) {
	double objective = 0;

	for (int i = 0; i < this.store.getShelves().size(); i++) {
	    Shelf shelf = this.store.getShelves().get(i);
	    int nh = shelf.getHorizontal();

	    for (int k = 0; k < shelf.getSegments().size(); k++) {
		Segment segment = shelf.getSegments().get(k);

		for (int j = 0; j < this.store.getProducts().size(); j++) {
		    Product product = this.store.getProducts().get(j);

		    double fsc = segment.getAttractiveness() * this.s[k][j] / segment.getCapacity();
		    objective += product.getMaxProfit() * fsc;
		    switch (obj) {
		    case AVA:
			objective -= param * this.y[k][j] / product.getHealthScore();
			break;
		    case HLUR:
			if (k % nh != 0) {
			    // Shelf segment k is not the rightmost segment.
			    double yhh = 0;

			    for (int k2 = k + 1; k2 < k + (nh - k % nh); k2++) {
				// Shelf segments to the right of shelf segment k.

				for (int j2 = 0; j2 < this.store.getProducts().size(); j2++) {
				    Product product2 = this.store.getProducts().get(j2);
				    double hh = product2.getHealthScore() - product.getHealthScore();
				    yhh += this.y[k2][j2] * hh;
				}
			    }
			    objective -= param * this.y[k][j] * yhh;
			}
			break;
		    default:
			objective -= param / product.getHealthScore() * fsc;
			break;
		    }
		}
	    }
	}
	return objective;
    }

    /**
     * Computes the objective value of this solution for the objective function of
     * type APSA.
     * 
     * @return the objective value
     */
    public double getObjectiveAPSA() {
	double objective = 0;

	for (int k = 0; k < this.store.getSegments().size(); k++) {
	    Segment segment = this.store.getSegments().get(k);

	    for (int j = 0; j < this.store.getProducts().size(); j++) {
		Product product = this.store.getProducts().get(j);

		double fsc = segment.getAttractiveness() * this.s[k][j] / segment.getCapacity();
		objective += product.getMaxProfit() * fsc;
	    }
	}
	return objective;
    }

    /**
     * Computes the objective value of this solution for the objective function of
     * type HAPSA.
     * 
     * @param gamma the visibility penalty parameter gamma
     * @param theta the healthy-left, unhealthy-right parameter theta
     * @return the objective value
     */
    public double getObjectiveHAPSA(double gamma, double theta) {
	double objective = 0;

	for (int i = 0; i < this.store.getShelves().size(); i++) {
	    Shelf shelf = this.store.getShelves().get(i);
	    int nh = shelf.getHorizontal();

	    for (int k = 0; k < shelf.getSegments().size(); k++) {
		Segment segment = shelf.getSegments().get(k);

		for (int j = 0; j < this.store.getProducts().size(); j++) {
		    Product product = this.store.getProducts().get(j);

		    double fsc = segment.getAttractiveness() * this.s[k][j] / segment.getCapacity();
		    objective += product.getMaxProfit() * fsc;

		    // Visibility penalty.
		    objective -= gamma / product.getHealthScore() * fsc;

		    // Healthy-left, unhealthy-right approach.
		    if (k % nh != 0) {
			// Shelf segment k is not the rightmost segment.
			double yhh = 0;

			for (int k2 = k + 1; k2 < k + (nh - k % nh); k2++) {
			    // Shelf segments to the right of shelf segment k.

			    for (int j2 = 0; j2 < this.store.getProducts().size(); j2++) {
				Product product2 = this.store.getProducts().get(j2);
				double hh = product2.getHealthScore() - product.getHealthScore();
				yhh += this.y[k2][j2] * hh;
			    }
			}
			objective -= theta * this.y[k][j] * yhh;
		    }
		}
	    }
	}
	return objective;
    }

    /**
     * Gets a partial solution, corresponding to a set of shelf indices.
     * 
     * @param shelves the array of shelf indices
     * @return the partial solution (resized)
     */
    public Solution getPartial(int[] shelves) {
	int ni = shelves.length;
	int nk = this.store.getShelves().get(0).getSegments().size();
	int nj = this.store.getProducts().size();
	double[][] newS = new double[ni * nk][nj];
	double[][] newX = new double[ni][nj];
	double[][] newY = new double[ni * nk][nj];
	ArrayList<Shelf> newShelves = new ArrayList<Shelf>(ni);

	for (int sh = 0; sh < ni; sh++) {
	    int i = shelves[sh];
	    newX[sh] = this.x[i].clone();
	    for (int k = 0; k < nk; k++) {
		newS[sh * nk + k] = this.s[i * nk + k].clone();
		newY[sh * nk + k] = this.y[i * nk + k].clone();
	    }
	    newShelves.add(this.store.getShelves().get(i));
	}

	ArrayList<Integer> considered = new ArrayList<Integer>(nj);
	for (int j = 0; j < nj; j++) {
	    considered.add(j);
	}
	Store newStore = this.partialStore(newShelves, considered);
	return new Solution(newS, newX, newY, newStore);
    }

    /**
     * Gets the decision variable s_kj.
     *
     * @return the decision variable s_kj
     */
    public double[][] getS() {
	return this.s;
    }

    public double getShelfObjective(Shelf shelf, Objective obj, double param) {
	int i = this.store.getShelves().indexOf(shelf);
	Solution partial = this.getPartial(new int[] { i });
	return partial.getObjective(obj, param);
    }

    public double getShelfObjectiveAPSA(Shelf shelf) {
	int i = this.store.getShelves().indexOf(shelf);
	Solution partial = this.getPartial(new int[] { i });
	return partial.getObjectiveAPSA();
    }

    public double getShelfObjectiveHAPSA(Shelf shelf, double gamma, double theta) {
	int i = this.store.getShelves().indexOf(shelf);
	Solution partial = this.getPartial(new int[] { i });
	return partial.getObjectiveHAPSA(gamma, theta);
    }

    /**
     * Gets the store.
     *
     * @return the store
     */
    public Store getStore() {
	return this.store;
    }

    /**
     * Gets the decision variable x_ij.
     *
     * @return the decision variable x_ij
     */
    public double[][] getX() {
	return this.x;
    }

    /**
     * Gets the decision variable y_kj.
     * 
     * @return the decision variable y_kj
     */
    public double[][] getY() {
	return this.y;
    }

    /**
     * Copies part of the store, based on a specified set of shelves. The products
     * in this partial store are either selected on this set of shelves or not
     * selected in the assortment at all.
     * 
     * @param shelves  the set of shelves to be included in the partial store
     * @param products the set of product indices that are to be considered
     * @return the partial copy of the store
     */
    public Store partialStore(ArrayList<Shelf> shelves, ArrayList<Integer> products) {
	int nSegments = shelves.size() * this.store.getShelves().get(0).getSegments().size();
	ArrayList<Segment> segments = new ArrayList<Segment>(nSegments);
	for (Shelf shelf : shelves) {
	    for (Segment segment : shelf.getSegments()) {
		segments.add(segment);
	    }
	}
	ArrayList<Product> partialProducts = new ArrayList<Product>();
	for (int j : products) {
	    partialProducts.add(this.store.getProducts().get(j));
	}

	// Remove all relationship pairs that contain non-considered products.
	HashSet<SymmetricPair> AllocationAffinity = new HashSet<SymmetricPair>();
	HashSet<SymmetricPair> AllocationDisaffinity = new HashSet<SymmetricPair>();
	HashSet<AsymmetricPair> AsymmetricAssortment = new HashSet<AsymmetricPair>();
	HashSet<SymmetricPair> SymmetricAssortment = new HashSet<SymmetricPair>();
	for (SymmetricPair pair : this.store.getAllocationAffinity()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (products.contains(j1) && products.contains(j2)) {
		int newj1 = products.indexOf(j1);
		int newj2 = products.indexOf(j2);
		SymmetricPair newPair = new SymmetricPair(newj1, pair.getProduct1(), newj2, pair.getProduct2());
		AllocationAffinity.add(newPair);
	    }
	}
	for (SymmetricPair pair : this.store.getAllocationDisaffinity()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (products.contains(j1) && products.contains(j2)) {
		int newj1 = products.indexOf(j1);
		int newj2 = products.indexOf(j2);
		SymmetricPair newPair = new SymmetricPair(newj1, pair.getProduct1(), newj2, pair.getProduct2());
		AllocationDisaffinity.add(newPair);
	    }
	}
	for (AsymmetricPair pair : this.store.getAsymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (products.contains(j1) && products.contains(j2)) {
		int newj1 = products.indexOf(j1);
		int newj2 = products.indexOf(j2);
		AsymmetricPair newPair = new AsymmetricPair(newj1, pair.getProduct1(), newj2, pair.getProduct2());
		AsymmetricAssortment.add(newPair);
	    }
	}
	for (SymmetricPair pair : this.store.getSymmetricAssortment()) {
	    int j1 = pair.getIndex1();
	    int j2 = pair.getIndex2();
	    if (products.contains(j1) && products.contains(j2)) {
		int newj1 = products.indexOf(j1);
		int newj2 = products.indexOf(j2);
		SymmetricPair newPair = new SymmetricPair(newj1, pair.getProduct1(), newj2, pair.getProduct2());
		SymmetricAssortment.add(newPair);
	    }
	}
	return new Store(partialProducts, shelves, segments, AllocationAffinity, AllocationDisaffinity,
		AsymmetricAssortment, SymmetricAssortment);
    }

    /**
     * Updates the current solution based on partial results, corresponding to a set
     * of shelf indices.
     * 
     * @param shelves  the array of shelf indices
     * @param products the list of indices of the considered products
     * @param newS     the new (partial) decision variable s_kj
     * @param newX     the new (partial) decision variable x_ij
     * @param newY     the new (partial) decision variable y_kj
     */
    public void updatePartial(int[] shelves, ArrayList<Integer> products, double[][] newS, double[][] newX,
	    double[][] newY) {
	int ni = shelves.length;
	int nk = this.store.getShelves().get(0).getSegments().size();
	int nj = products.size();

	for (int sh = 0; sh < ni; sh++) {
	    int i = shelves[sh];
	    for (int pr = 0; pr < nj; pr++) {
		int j = products.get(pr);
		this.x[i][j] = newX[sh][pr];
		for (int k = 0; k < nk; k++) {
		    this.s[i * nk + k][j] = newS[sh * nk + k][pr];
		    this.y[i * nk + k][j] = newS[sh * nk + k][pr];
		}
	    }
	}
    }
}
