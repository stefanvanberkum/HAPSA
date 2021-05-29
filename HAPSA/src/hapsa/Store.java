/*
 * Store.java
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
import java.util.HashSet;

/**
 * The Store class provides definitions and methods for a store.
 *
 * @author Stefan van Berkum
 *
 */
public class Store {

    /** The set of product category pairs that have allocation affinity. */
    private HashSet<SymmetricPair> allocationAffinity;

    /** The set of product category pairs that have allocation disaffinity. */
    private HashSet<SymmetricPair> allocationDisaffinity;

    /**
     * The set of product category pairs that have asymmetric assortment affinity.
     */
    private HashSet<AsymmetricPair> asymmetricAssortment;

    /** The set of product categories */
    private ArrayList<Product> products;

    /** The set of shelf segments in this store. */
    private ArrayList<Segment> segments;

    /** The set of shelves in this store. */
    private ArrayList<Shelf> shelves;

    /**
     * The set of product category pairs that have symmetric assortment affinity.
     */
    private HashSet<SymmetricPair> symmetricAssortment;

    /**
     * Constructs a store.
     * 
     * @param products              the set of products considered in this store
     * @param shelves               the set of shelves in this store
     * @param segments              the set segments in this store
     * @param allocationAffinity    the set of product pairs that have allocation
     *                              affinity
     * @param allocationDisaffinity the set of product pairs that have allocation
     *                              disaffinity
     * @param asymmetricAssortment  the set of product pairs that have asymmetric
     *                              assortment affinity
     * @param symmetricAssortment   set of product pairs that have symmetric
     *                              assortment affinity
     */
    public Store(ArrayList<Product> products, ArrayList<Shelf> shelves, ArrayList<Segment> segments,
	    HashSet<SymmetricPair> allocationAffinity, HashSet<SymmetricPair> allocationDisaffinity,
	    HashSet<AsymmetricPair> asymmetricAssortment, HashSet<SymmetricPair> symmetricAssortment) {
	this.products = products;
	this.shelves = shelves;
	this.segments = segments;
	this.allocationAffinity = allocationAffinity;
	this.allocationDisaffinity = allocationDisaffinity;
	this.asymmetricAssortment = asymmetricAssortment;
	this.symmetricAssortment = symmetricAssortment;
    }

    /**
     * Gets the set of product pairs that have allocation affinity.
     *
     * @return the set of product pairs that have allocation affinity
     */
    public HashSet<SymmetricPair> getAllocationAffinity() {
	return this.allocationAffinity;
    }

    /**
     * Gets the set of product pairs that have allocation disaffinity.
     *
     * @return the set of product pairs that have allocation disaffinity
     */
    public HashSet<SymmetricPair> getAllocationDisaffinity() {
	return this.allocationDisaffinity;
    }

    /**
     * Gets the set of product pairs that have asymmetric assortment affinity.
     *
     * @return the set of product pairs that have asymmetric assortment affinity
     */
    public HashSet<AsymmetricPair> getAsymmetricAssortment() {
	return this.asymmetricAssortment;
    }

    /**
     * Gets the products.
     *
     * @return the products
     */
    public ArrayList<Product> getProducts() {
	return this.products;
    }

    /**
     * Gets the segments.
     *
     * @return the segments
     */
    public ArrayList<Segment> getSegments() {
	return this.segments;
    }

    /**
     * Gets the shelves.
     *
     * @return the shelves
     */
    public ArrayList<Shelf> getShelves() {
	return this.shelves;
    }

    /**
     * Gets the set of product pairs that have symmetric assortment affinity.
     *
     * @return the set of product pairs that have symmetric assortment affinity
     */
    public HashSet<SymmetricPair> getSymmetricAssortment() {
	return this.symmetricAssortment;
    }
}
