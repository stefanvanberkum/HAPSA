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

import java.io.File;

import ilog.concert.IloIntVar;
import ilog.concert.IloNumVar;

/**
 * The Utils class provides some general utility methods for the HAPSA method.
 *
 * @author Stefan van Berkum
 *
 */
public class Utils {

    /**
     * Concatenates a matrix of decision variables into a vector.
     * 
     * @param var the matrix to be concatenated
     * @return the concatenated vector
     */
    static IloNumVar[] concatenate(IloNumVar[][] var) {
	IloNumVar[] result = new IloNumVar[var.length * var[0].length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		result[i * var[0].length + j] = var[i][j];
	    }
	}
	return result;
    }

    /**
     * Concatenates two matrices of decision variables into a single vector.
     * 
     * @param var1 the first matrix to be concatenated
     * @param var2 the second matrix to be concatenated
     * @return the concatenated vector
     */
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

    /**
     * Concatenates three matrices of decision variables into a single vector.
     * 
     * @param var1 the first matrix to be concatenated
     * @param var2 the second matrix to be concatenated
     * @param var3 the third matrix to be concatenated
     * @return the concatenated vector
     */
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

    /**
     * Gets the transpose of a matrix of decision variables.
     * 
     * @param var the matrix to be transposed
     * @return the transposed matrix
     */
    static IloIntVar[][] getTranspose(IloIntVar[][] var) {
	IloIntVar[][] transpose = new IloIntVar[var[0].length][var.length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		transpose[j][i] = var[i][j];
	    }
	}
	return transpose;
    }

    /**
     * Gets the transpose of a matrix of decision variables.
     * 
     * @param var the matrix to be transposed
     * @return the transposed matrix
     */
    static IloNumVar[][] getTranspose(IloNumVar[][] var) {
	IloNumVar[][] transpose = new IloNumVar[var[0].length][var.length];
	for (int i = 0; i < var.length; i++) {
	    for (int j = 0; j < var[0].length; j++) {
		transpose[j][i] = var[i][j];
	    }
	}
	return transpose;
    }

    /**
     * Assigns a version number to this run for objective function types: AVA, HLUR,
     * and VIS.
     * 
     * @param obj      the objective function type, one of: AVA, HLUR, or VIS
     * @param shelves  the number of shelves used for this run
     * @param products the number of product used for this run
     * @param param    the parameter associated with the objective function type
     * @param version  the version number
     */
    static void setVersion(Model.Objective obj, int shelves, int products, double param, int version) {
	String path = "Results/" + obj + "/";

	String pathCSV = path + "/CSV/" + shelves + "_" + products + "_" + param;
	File oldScoresCSV = new File(pathCSV + "_scores.csv");
	File newScoresCSV = new File(pathCSV + "_scores_v" + version + ".csv");
	oldScoresCSV.renameTo(newScoresCSV);
	File oldShelfCSV = new File(pathCSV + "_shelf.csv");
	File newShelfCSV = new File(pathCSV + "_shelf_v" + version + ".csv");
	oldShelfCSV.renameTo(newShelfCSV);

	String pathSummary = path + "/Summary/" + shelves + "_" + products + "_" + param;
	File oldSummary = new File(pathSummary + ".txt");
	File newSummary = new File(pathSummary + "_v" + version + ".txt");
	oldSummary.renameTo(newSummary);

	String[] vars = new String[] { "s", "x", "y" };
	for (String var : vars) {
	    String pathVar = path + "/Variables/" + shelves + "_" + products + "_" + param + "_" + var;
	    File oldVar = new File(pathVar + ".csv");
	    File newVar = new File(pathVar + "_v" + version + ".csv");
	    oldVar.renameTo(newVar);
	}
    }

    /**
     * Assigns a version number to this run for objective function type APSA.
     * 
     * @param shelves  the number of shelves used for this run
     * @param products the number of product used for this run
     * @param version  the version number
     */
    static void setVersionAPSA(int shelves, int products, int version) {
	String path = "Results/HAPSA/";

	String pathCSV = path + "/CSV/" + shelves + "_" + products;
	File oldScoresCSV = new File(pathCSV + "_scores.csv");
	File newScoresCSV = new File(pathCSV + "_scores_v" + version + ".csv");
	oldScoresCSV.renameTo(newScoresCSV);
	File oldShelfCSV = new File(pathCSV + "_shelf.csv");
	File newShelfCSV = new File(pathCSV + "_shelf_v" + version + ".csv");
	oldShelfCSV.renameTo(newShelfCSV);

	String pathSummary = path + "/Summary/" + shelves + "_" + products;
	File oldSummary = new File(pathSummary + ".txt");
	File newSummary = new File(pathSummary + "_v" + version + ".txt");
	oldSummary.renameTo(newSummary);

	String[] vars = new String[] { "s", "x", "y" };
	for (String var : vars) {
	    String pathVar = path + "/Variables/" + shelves + "_" + products + "_" + var;
	    File oldVar = new File(pathVar + ".csv");
	    File newVar = new File(pathVar + "_v" + version + ".csv");
	    oldVar.renameTo(newVar);
	}
    }

    /**
     * Assigns a version number to this run for objective function type HAPSA.
     * 
     * @param shelves  the number of shelves used for this run
     * @param products the number of product used for this run
     * @param gamma    the gamma parameter for the visibility penalty
     * @param theta    the theta parameter for the healthy-left, unhealthy-right
     *                 approach
     * @param version  the version number
     */
    static void setVersionHAPSA(int shelves, int products, double gamma, double theta, int version) {
	String path = "Results/HAPSA/";

	String pathCSV = path + "/CSV/" + shelves + "_" + products + "_" + gamma + "_" + theta;
	File oldScoresCSV = new File(pathCSV + "_scores.csv");
	File newScoresCSV = new File(pathCSV + "_scores_v" + version + ".csv");
	oldScoresCSV.renameTo(newScoresCSV);
	File oldShelfCSV = new File(pathCSV + "_shelf.csv");
	File newShelfCSV = new File(pathCSV + "_shelf_v" + version + ".csv");
	oldShelfCSV.renameTo(newShelfCSV);

	String pathSummary = path + "/Summary/" + shelves + "_" + products + "_" + gamma + "_" + theta;
	File oldSummary = new File(pathSummary + ".txt");
	File newSummary = new File(pathSummary + "_v" + version + ".txt");
	oldSummary.renameTo(newSummary);

	String[] vars = new String[] { "s", "x", "y" };
	for (String var : vars) {
	    String pathVar = path + "/Variables/" + shelves + "_" + products + "_" + gamma + "_" + theta + "_" + var;
	    File oldVar = new File(pathVar + ".csv");
	    File newVar = new File(pathVar + "_v" + version + ".csv");
	    oldVar.renameTo(newVar);
	}
    }
}
