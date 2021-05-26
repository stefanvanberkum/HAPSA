/*
 * Product.java
 * 
 * Version 1.0
 * 
 * 12/05/2021
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

/**
 * The Product class provides definitions and methods for a product within the
 * context of a store.
 * 
 * @author Stefan van Berkum
 *
 */
public class Product {

    /** The expected demand of this product on a single segment. */
    private int demand;

    /** True if this product is a fast mover. */
    private boolean fastMover;

    /** The health score of this product. */
    private int healthScore;

    /**
     * The impulse purchase potential of this product, between zero (exclusive) and
     * one (inclusive).
     */
    private double impulse;

    /** The profit margin of this product. */
    private double margin;

    /**
     * The maximum attainable profit on a single segment for this product. Equals
     * {@link #margin} * {@link #impulse} * {@link #demand}.
     */
    private double maxProfit;

    /** The maximum space to be allocated to this product. */
    private double maxSpace;

    /**
     * The minimum space to be allocated on any segment that this product is
     * assigned to.
     */
    private double minAllocated;

    /** The minimum space to be allocated to this product. */
    private double minSpace;

    /**
     * Constructs a product. Computes maximum attainable profit with given profit
     * margin, impulse purchase potential, and expected demand.
     * 
     * @param fastMover    True if this product is a fast mover
     * @param margin       The profit margin of this product
     * @param impulse      The impulse purchase potential of this product, between
     *                     zero (exclusive) and one (inclusive)
     * @param minSpace     The minimum space to be allocated to this product
     * @param maxSpace     The maximum space to be allocated to this product
     * @param minAllocated The minimum space to be allocated on any segment that
     *                     this product is assigned to
     * @param demand       The expected demand of this product on a single segment
     */
    public Product(boolean fastMover, double margin, double impulse, double minSpace, double maxSpace,
	    double minAllocated, int demand) {
	this.fastMover = fastMover;
	this.margin = margin;
	this.impulse = impulse;
	this.minSpace = minSpace;
	this.maxSpace = maxSpace;
	this.minAllocated = minAllocated;
	this.demand = demand;
	this.maxProfit = margin * impulse * demand;
    }

    /**
     * Constructs a product. Takes manual maximum attainable profit.
     * 
     * @param fastMover    True if this product is a fast mover
     * @param maxProfit    The maximum attainable profit on a single segment for
     *                     this product
     * @param minSpace     The minimum space to be allocated to this product
     * @param maxSpace     The maximum space to be allocated to this product
     * @param minAllocated The minimum space to be allocated on any segment that
     *                     this product is assigned to
     * @param healthScore  The health score of this product
     */
    public Product(boolean fastMover, double maxProfit, double minSpace, double maxSpace, double minAllocated,
	    int healthScore) {
	this.fastMover = fastMover;
	this.maxProfit = maxProfit;
	this.minSpace = minSpace;
	this.maxSpace = maxSpace;
	this.minAllocated = minAllocated;
	this.healthScore = healthScore;
    }

    /**
     * Gets the demand.
     *
     * @return the demand
     */
    public int getDemand() {
	return this.demand;
    }

    /**
     * Gets the health score.
     *
     * @return the health score
     */
    public int getHealthScore() {
	return this.healthScore;
    }

    /**
     * Gets the impulse.
     *
     * @return the impulse
     */
    public double getImpulse() {
	return this.impulse;
    }

    /**
     * Gets the margin.
     *
     * @return the margin
     */
    public double getMargin() {
	return this.margin;
    }

    /**
     * Gets the maxProfit.
     *
     * @return the maxProfit
     */
    public double getMaxProfit() {
	return this.maxProfit;
    }

    /**
     * Gets the maxSpace.
     *
     * @return the maxSpace
     */
    public double getMaxSpace() {
	return this.maxSpace;
    }

    /**
     * Gets the minAllocated.
     *
     * @return the minAllocated
     */
    public double getMinAllocated() {
	return this.minAllocated;
    }

    /**
     * Gets the minSpace.
     *
     * @return the minSpace
     */
    public double getMinSpace() {
	return this.minSpace;
    }

    /**
     * Gets the fastMover.
     *
     * @return the fastMover
     */
    public boolean isFastMover() {
	return this.fastMover;
    }
}
