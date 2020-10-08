package Sudoku;

import java.util.*;

public class Sudoku {
	private int[][] matrix;

	/** Class constructor.
	 * 
	 * Creates Sudoku with 9x9 board. */
	public Sudoku() {
		matrix = new int[9][9];
	}

	/** Adds number to row r, column c. 
	 * 
	 * @param r row
	 * @param c column
	 * @param number value
	 * */
	public void add(int r, int c, int number) {
			matrix[r][c] = number;
	}

	/** Returns the value at row r, column c.
	 * 
	 *  @param r row
	 *  @param c column
	 *  @return value at r,c
	 *  */
	public int get(int r, int c) {
		return matrix[r][c];
	}

	/** Clears the Sudoku grid.*/
	public void clear() {
		for (int r = 0; r < 9; r++) {
			for (int c = 0; c < 9; c++) {
				matrix[r][c] = 0;
			}
		}
	}

	/** Checks if the user's board is ok.
	 *  
	 *  @return true if user's board is valid.
	 * 
	 */
	public boolean checkBoard() {
		for(int i = 0; i<9; i ++) {
			for(int j = 0; j<9; j++) {
				if (matrix[i][j] != 0) {
					int n = matrix[i][j];
					matrix[i][j] = 0;
					boolean b1 = checkRow(i, n);
					boolean b2 = checkColumn(j, n);
					boolean b3 = checkBox(i, j, n);
					matrix[i][j] = n;
					if(b1 == false || b2 == false || b3 == false) {
						return false;
					}
				}
			}
		}
		return true;
	}
	
	
	/** Solves the sudoku with backtracking
	 *  
	 * @return true if the sudoku is solvable
	 * 
	 * */
	public boolean solve() {
		if(checkBoard()) {
		return solve(0, 0);
		} else {
			return false;
		}
	}

	/** Solves the sudoku with backtracking */
	private boolean solve(int r, int c) {

		if (matrix[r][c] == 0) {

			for (int i = 1; i < 10; i++) {
				boolean b1 = checkRow(r, i);
				boolean b2 = checkColumn(c, i);
				boolean b3 = checkBox(r, c, i);
				if (b1 && b2 && b3) {
					add(r, c, i);

					// move to next square
					int ctemp = c + 1;
					int rtemp = r;
					if (ctemp == 9) {
						ctemp = 0;
						rtemp++;
					}
					if (rtemp == 9) {
						return true;
					}

					if (solve(rtemp, ctemp)) {
						return true;
					}
				}
			}
			matrix[r][c] = 0;
			return false;
		} else {
			int ctemp = c + 1;
			int rtemp = r;
			if (ctemp == 9) {
				ctemp = 0;
				rtemp++;
			}
			if (rtemp == 9) {
				return true;
			}
			return solve(rtemp, ctemp);
		}
	}

	/** Checks if there is a square with the number nbr in the row r 
	 * 
	 * @param r row
	 * @param nbr value
	 * @return true if the value can be put on the row
	 * */
	private boolean checkRow(int r, int nbr) { // nbr än inte inlaggt i matrix än.
		for (int col = 0; col < 9; col++) {
			if (nbr == matrix[r][col]) {
				return false;
			}
		}
		return true;
	}

	/** Checks if there is a square with the number nbr in the column c 
	 * 
	 * @param c column
	 * @param nbr value
	 * @return true if the value can be put in the column
	 * 
	 * */
	private boolean checkColumn(int c, int nbr) { // nbr än inte inlaggt i matrix än.
		for (int row = 0; row < 9; row++) {
			if (nbr == matrix[row][c]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Checks if there is a square with the number nbr in the box that contains the
	 * square with row r, column c
	 * 
	 * @param r row
	 * @param c column
	 * @param nbr value
	 * @return true if the value can be put in the box
	 * 
	 */
	private boolean checkBox(int r, int c, int nbr) { // nbr än inte inlaggt i matrix än.
		if (c < 3) { // box 1, 4, 7.
			if (r < 3) { // box 1
				for (int row = 0; row < 3; row++) {
					for (int col = 0; col < 3; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 6) { // box 4
				for (int row = 3; row < 6; row++) {
					for (int col = 0; col < 3; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 9) { // box 7
				for (int row = 6; row < 9; row++) {
					for (int col = 0; col < 3; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			}
		} else if (c < 6) { // box 2, 5, 8.
			if (r < 3) { // box
				for (int row = 0; row < 3; row++) {
					for (int col = 3; col < 6; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 6) { // box 5
				for (int row = 3; row < 6; row++) {
					for (int col = 3; col < 6; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 9) { // box 8
				for (int row = 6; row < 9; row++) {
					for (int col = 3; col < 6; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			}
		} else if (c < 9) { // box 3, 6, 9.
			if (r < 3) { // box 3
				for (int row = 0; row < 3; row++) {
					for (int col = 6; col < 9; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 6) { // box 6
				for (int row = 3; row < 6; row++) {
					for (int col = 6; col < 9; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			} else if (r < 9) { // box 9
				for (int row = 6; row < 9; row++) {
					for (int col = 6; col < 9; col++) {
						if (nbr == matrix[row][col]) {
							return false;
						}
					}
				}
			}
		}
		return true;
	}
}
