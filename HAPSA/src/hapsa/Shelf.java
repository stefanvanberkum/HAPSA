/*
 * Shelf.java
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

import java.util.ArrayList;

/**
 * The Shelf class provides definitions and methods for a store shelf. Shelves
 * are sorted based on their relative attractiveness (non-increasing order).
 *
 * @author Stefan van Berkum
 *
 */
public class Shelf implements Comparable<Shelf> {

    /** The relative attractiveness of this shelf. */
    private double attractiveness = 0;

    /** The number of horizontal positions on this shelf. */
    private int horizontal;

    /** The set of segments on this shelf. */
    private ArrayList<Segment> segments;

    /**
     * Constructs a store shelf.
     * 
     * @param segments   the set of segments on this shelf
     * @param horizontal the number of horizontal positions on this shelf
     */
    public Shelf(ArrayList<Segment> segments, int horizontal) {
	this.segments = segments;
	this.horizontal = horizontal;

	double totalAttr = 0;
	double totalCap = 0;
	for (Segment segment : segments) {
	    totalAttr += segment.getAttractiveness() * segment.getCapacity();
	    totalCap += segment.getCapacity();
	}
	this.attractiveness = totalAttr / totalCap;
    }

    @Override
    public int compareTo(Shelf otherShelf) {
	return Double.compare(otherShelf.attractiveness, this.attractiveness);
    }

    /**
     * Gets the number of horizontal positions.
     *
     * @return the number horizontal positions
     */
    public int getHorizontal() {
	return this.horizontal;
    }

    /**
     * Gets the segments.
     *
     * @return the segments
     */
    public ArrayList<Segment> getSegments() {
	return this.segments;
    }
}
