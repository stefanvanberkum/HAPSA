/*
 * Utils.java
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

import ilog.concert.IloIntVar;
import ilog.concert.IloNumVar;

/**
 * The Utils class provides some general utility methods for the HAPSA method.
 *
 * @author Stefan van Berkum
 *
 */
public class Utils {

    static IloNumVar[] concatenate(IloNumVar[][] var) {
	IloNumVar[] result = new IloNumVar[var.length * var[0].length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		result[i * var[0].length + j] = var[i][j];
	    }
	}
	return result;
    }

    static IloNumVar[] concatenate(IloNumVar[][] var1, IloNumVar[][] var2) {
	IloNumVar[] concat1 = concatenate(var1);
	IloNumVar[] concat2 = concatenate(var2);
	IloNumVar[] result = new IloNumVar[concat1.length + concat2.length];
	int i = 0;
	for (int j = 0; j < concat1.length; j++) {
	    result[i++] = concat1[j];
	}
	for (int j = 0; j < concat2.length; j++) {
	    result[i++] = concat2[j];
	}
	return result;
    }

    static IloNumVar[] concatenate(IloNumVar[][] var1, IloNumVar[][] var2, IloNumVar[][] var3) {
	IloNumVar[] oneTwo = concatenate(var1, var2);
	IloNumVar[] three = concatenate(var3);
	IloNumVar[] result = new IloNumVar[oneTwo.length + three.length];
	int i = 0;
	for (int j = 0; j < oneTwo.length; j++) {
	    result[i++] = oneTwo[j];
	}
	for (int j = 0; j < three.length; j++) {
	    result[i++] = three[j];
	}
	return result;
    }

    static IloIntVar[][] getTranspose(IloIntVar[][] var) {
	IloIntVar[][] transpose = new IloIntVar[var[0].length][var.length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		transpose[j][i] = var[i][j];
	    }
	}
	return transpose;
    }

    static IloNumVar[][] getTranspose(IloNumVar[][] var) {
	IloNumVar[][] transpose = new IloNumVar[var[0].length][var.length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		transpose[j][i] = var[i][j];
	    }
	}
	return transpose;
    }
}
