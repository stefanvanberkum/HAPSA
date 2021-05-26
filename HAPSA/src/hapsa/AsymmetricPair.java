/*
 * AsymmetricPair.java
 * 
 * v1.0
 *
 * 15 May 2021
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
 * The AsymmetricPair class defines a basic framework for a pair of products
 * which have an asymmetric affinity relationship.
 *
 * @author Stefan van Berkum
 *
 */
public class AsymmetricPair {

    /** The index j of the first product. */
    private int index1;

    /** The index j' of the second product. */
    private int index2;

    /** The first product. */
    private Product product1;

    /** The second product. */
    private Product product2;

    /**
     * Constructs a product pair with an asymmetric affinity relationship.
     * 
     * @param index1   the index of the first product
     * @param product1 the first product
     * @param index2   the index of the second product
     * @param product2 the second product
     * @throws IllegalArgumentException if one of the passed products is null
     */
    public AsymmetricPair(int index1, Product product1, int index2, Product product2) throws IllegalArgumentException {
	if (product1 == null || product2 == null) {
	    throw new IllegalArgumentException("Invalid pair, both products must be non-null.");
	}
	this.index1 = index1;
	this.product1 = product1;
	this.index2 = index2;
	this.product2 = product2;
    }

    @Override
    public boolean equals(Object obj) {
	if (this == obj) {
	    return true;
	}
	if (!(obj instanceof AsymmetricPair)) {
	    return false;
	}

	AsymmetricPair other = (AsymmetricPair) obj;
	if (this.product1 == null || this.product2 == null || other.product1 == null || other.product2 == null) {
	    return false;
	}
	return (this.product1.equals(other.product1) && this.product2.equals(other.product2));
    }

    /**
     * Gets the first index.
     *
     * @return the fist index
     */
    public int getIndex1() {
	return this.index1;
    }

    /**
     * Gets the second index.
     *
     * @return the second index
     */
    public int getIndex2() {
	return this.index2;
    }

    /**
     * Gets the first product.
     *
     * @return the first product
     */
    public Product getProduct1() {
	return this.product1;
    }

    /**
     * Gets the second product.
     *
     * @return the second product
     */
    public Product getProduct2() {
	return this.product2;
    }

    @Override
    public int hashCode() {
	final int prime = 31;
	int result = 1;
	result = prime * result + ((this.product1 == null) ? 0 : this.product1.hashCode());
	result = prime * result + ((this.product2 == null) ? 0 : this.product2.hashCode());
	return result;
    }
}
