/*
 * Segment.java
 * 
 * v1.0
 *
 * 12 May 2021
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
 * The Segment class provides definitions and methods for a shelf segment.
 * Segments are sorted based on their capacity.
 *
 * @author Stefan van Berkum
 *
 */
public class Segment implements Comparable<Segment> {

    /**
     * The attractiveness score of this segment, between zero (exclusive) and one
     * (inclusive).
     */
    private double attractiveness;

    /** The capacity of this segment. */
    private double capacity;

    /**
     * Constructs a shelf segment.
     * 
     * @param attractiveness the attractiveness score of this segment, between zero
     *                       (exclusive) and one (inclusive)
     * @param capacity       the capacity of this segment
     */
    public Segment(double attractiveness, double capacity) {
	this.attractiveness = attractiveness;
	this.capacity = capacity;
    }

    @Override
    public int compareTo(Segment otherSegment) {
	return Double.compare(this.capacity, otherSegment.capacity);
    }

    /**
     * Gets the attractiveness.
     *
     * @return the attractiveness
     */
    public double getAttractiveness() {
	return this.attractiveness;
    }

    /**
     * Gets the capacity.
     *
     * @return the capacity
     */
    public double getCapacity() {
	return this.capacity;
    }
}
