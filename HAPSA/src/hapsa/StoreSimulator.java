/*
 * StoreSimulator.java
 * 
 * v1.0
 *
 * 17 May 2021
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
import java.util.Random;

/**
 * The StoreSimulator class provides methods for the simulation of a store.
 *
 * @author Stefan van Berkum
 *
 */
public class StoreSimulator {

    /** The capacity of each shelf segment. */
    private static double CAPACITY = 1;

    /** The horizontal attractiveness bonus for end-of-aisle segments. */
    private static double[] END_SEG = new double[] { 0.06, 0.1 };

    /** The lower bound for the randomly generated health score for each product. */
    private static double HEALTH_SCORE_LB = 1;

    /** The upper bound for the randomly generated health score for each product. */
    private static double HEALTH_SCORE_UB = 100;

    /** The number of horizontal positions on each shelf. */
    private static int HORIZONTAL = 3;

    /** The horizontal base shelf attractiveness categories (t). */
    private static double[] HORIZONTAL_CAT = new double[] { 0.05, 0.25, 0.45, 0.65, 0.85 };

    /**
     * The lower bound for the randomly generated maximum attainable profit for each
     * product.
     */
    private static double MAX_PROFIT_LB = 1;

    /**
     * The upper bound for the randomly generated maximum attainable profit for each
     * product.
     */
    private static double MAX_PROFIT_UB = 25;

    /**
     * The upper bound for the randomly generated maximum space for each product.
     */
    private static double MAX_SPACE_UB = 6;

    /** The horizontal attractiveness bonus for middle segments. */
    private static double[] MIDDLE_SEG = new double[] { 0.0, 0.05 };

    /** The minimum allocated space for each product (if selected). */
    private static double MIN_ALLOCATED = 0.1;

    /**
     * The lower bound for the randomly generated minimum space for each product.
     */
    private static double MIN_SPACE_LB = 1;

    /**
     * The upper bound for the randomly generated minimum space for each product.
     */
    private static double MIN_SPACE_UB = 3;

    /**
     * The fraction of products to be randomly selected for each affinity
     * relationship.
     */
    private static double RELATION_FRAC = 0.01;

    /** The number of vertical positions on each shelf. */
    private static int VERTICAL = 6;

    /** The vertical attractiveness categories, from the top to the bottom shelf. */
    private static double[] VERTICAL_CAT = new double[] { 0.925, 0.95, 1.0, 1.0, 0.95, 0.9 };

    /** The number of product categories. */
    private int nProducts;

    /** The number of shelves. */
    private int nShelves;

    /** The random number generator. */
    private Random rand;

    public StoreSimulator(int nProducts, int nShelves, int seed) {
	this.nProducts = nProducts;
	this.nShelves = nShelves;
	this.rand = new Random(seed);
    }

    /**
     * Simulates a store.
     * 
     * @return the simulated store
     */
    public Store simulate() {
	ArrayList<Product> products = new ArrayList<Product>(this.nProducts);
	ArrayList<Shelf> shelves = new ArrayList<Shelf>(this.nShelves);
	ArrayList<Segment> segments = new ArrayList<Segment>(this.nShelves * HORIZONTAL * VERTICAL);

	int nSelect = Math.toIntExact(Math.round(RELATION_FRAC * this.nProducts));
	HashSet<SymmetricPair> allocationAffinity = new HashSet<SymmetricPair>(nSelect);
	HashSet<SymmetricPair> allocationDisaffinity = new HashSet<SymmetricPair>(nSelect);
	HashSet<AsymmetricPair> asymmetricAssortment = new HashSet<AsymmetricPair>(nSelect);
	HashSet<SymmetricPair> symmetricAssortment = new HashSet<SymmetricPair>(nSelect);

	// Simulate the product categories.
	for (int j = 0; j < this.nProducts; j++) {
	    products.add(this.simulateProduct());
	}

	// Simulate the shelves and segments.
	for (int i = 0; i < this.nShelves; i++) {
	    double base_attr = HORIZONTAL_CAT[i % HORIZONTAL_CAT.length];
	    Shelf shelf = this.simulateShelf(base_attr);
	    shelves.add(shelf);

	    for (Segment segment : shelf.getSegments()) {
		segments.add(segment);
	    }
	}

	// Simulate affinity relationships.
	for (int x = 0; x < nSelect; x++) {
	    allocationAffinity.add(this.randomSymmetric(products));
	    allocationDisaffinity.add(this.randomSymmetric(products));
	    asymmetricAssortment.add(this.randomAsymmetric(products));
	    symmetricAssortment.add(this.randomSymmetric(products));
	}
	return new Store(products, shelves, segments, allocationAffinity, allocationDisaffinity, asymmetricAssortment,
		symmetricAssortment);
    }

    /**
     * Generates a shelf segment with a given horizontal attractiveness score.
     * 
     * @param h_attr the horizontal attractiveness score of this shelf segment
     *               (between zero and one)
     * @param v_pos  the vertical position of this shelf segment (between zero and
     *               VERTICAL, where zero is the top shelf)
     * @return the generated shelf segment
     */
    private Segment generateSegment(double h_attr, int v_pos) {
	double v_attr = VERTICAL_CAT[v_pos];
	double attractiveness = h_attr * v_attr;
	return new Segment(attractiveness, CAPACITY);
    }

    /**
     * Generates a random asymmetric pair from a given list of products.
     * 
     * @param products the list of products
     * @return the random asymmetric pair
     */
    private AsymmetricPair randomAsymmetric(ArrayList<Product> products) {
	int index1 = this.rand.nextInt(products.size());
	Product product1 = products.get(index1);
	int index2 = this.rand.nextInt(products.size());
	Product product2 = products.get(index2);

	// Make sure the second product is not equal to the first product.
	while (product1 == product2) {
	    index2 = this.rand.nextInt(products.size());
	    product2 = products.get(index2);
	}

	try {
	    return new AsymmetricPair(index1, product1, index2, product2);
	} catch (IllegalArgumentException e) {
	    e.printStackTrace();
	    System.exit(-1);
	}
	return null;
    }

    /**
     * Generates a random symmetric pair from a given list of products.
     * 
     * @param products the list of products
     * @return the random symmetric pair
     */
    private SymmetricPair randomSymmetric(ArrayList<Product> products) {
	int index1 = this.rand.nextInt(products.size());
	Product product1 = products.get(index1);
	int index2 = this.rand.nextInt(products.size());
	Product product2 = products.get(index2);

	// Make sure the second product is not equal to the first product.
	while (product1 == product2) {
	    index2 = this.rand.nextInt(products.size());
	    product2 = products.get(index2);
	}

	try {
	    return new SymmetricPair(index1, product1, index2, product2);
	} catch (IllegalArgumentException e) {
	    e.printStackTrace();
	    System.exit(-1);
	}
	return null;
    }

    /**
     * Simulates a product category.
     * 
     * @return the simulated product category
     */
    private Product simulateProduct() {
	double randMaxProfit = MAX_PROFIT_LB + (MAX_PROFIT_UB - MAX_PROFIT_LB) * this.rand.nextDouble();
	double maxProfit = Math.round(100.0 * randMaxProfit) / 100.0;
	double randMinSpace = MIN_SPACE_LB + (MIN_SPACE_UB - MIN_SPACE_LB) * this.rand.nextDouble();
	int minSpace = Math.toIntExact(Math.round(randMinSpace));
	int maxSpace = minSpace;
	while (maxSpace == minSpace) {
	    double randMaxSpace = minSpace + (MAX_SPACE_UB - minSpace) * this.rand.nextDouble();
	    maxSpace = Math.toIntExact(Math.round(randMaxSpace));
	}
	double randHealthScore = HEALTH_SCORE_LB + (HEALTH_SCORE_UB - HEALTH_SCORE_LB) * this.rand.nextDouble();
	int healthScore = Math.toIntExact(Math.round(randHealthScore));
	return new Product(false, maxProfit, minSpace, maxSpace, MIN_ALLOCATED, healthScore);
    }

    /**
     * Simulates a shelf segment for a shelf with a given horizontal base
     * attractiveness score.
     * 
     * @param base_attr the horizontal base attractiveness score of the shelf
     *                  (between zero and one)
     * @param h_pos     the horizontal position of the segment (between zero and
     *                  HORIZONTAL)
     * @param v_pos     the vertical position of the segment (between zero and
     *                  VERTICAL, where zero is the top shelf)
     * @return the simulated segment
     */
    private Segment simulateSegment(double base_attr, int h_pos, int v_pos) {
	double h_attr;
	if (h_pos == 0 || h_pos == HORIZONTAL) {
	    // End-of-aisle shelf segment.
	    h_attr = base_attr + (END_SEG[0] + (END_SEG[1] - END_SEG[0]) * this.rand.nextDouble());
	} else {
	    // Middle shelf segment.
	    h_attr = base_attr + (MIDDLE_SEG[0] + (MIDDLE_SEG[1] - MIDDLE_SEG[0]) * this.rand.nextDouble());
	}
	double v_attr = VERTICAL_CAT[v_pos];
	double attractiveness = h_attr * v_attr;
	return new Segment(attractiveness, CAPACITY);
    }

    /**
     * Simulates a shelf with a given horizontal base attractiveness score.
     * 
     * @param base_attr the horizontal base attractiveness score
     * @return the simulated shelf
     */
    private Shelf simulateShelf(double base_attr) {
	ArrayList<Segment> segments = new ArrayList<Segment>(HORIZONTAL * VERTICAL);

	// Simulate upper shelf segments.
	for (int k = 0; k < HORIZONTAL; k++) {
	    segments.add(this.simulateSegment(base_attr, k, 0));
	}

	// Generate remaining shelf segments.
	for (int row = 1; row < VERTICAL; row++) {
	    for (int col = 0; col < HORIZONTAL; col++) {
		double h_attr = segments.get(col % HORIZONTAL).getAttractiveness() / VERTICAL_CAT[0];
		segments.add(this.generateSegment(h_attr, row));
	    }
	}
	return new Shelf(segments, HORIZONTAL);
    }
}
