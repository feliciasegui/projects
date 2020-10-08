package Sudoku;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import Sudoku.Sudoku;

public class TestSudoku {
	private Sudoku sudoku;

	@Before
	public void setUp() throws Exception {
		sudoku = new Sudoku();
	}

	@After
	public void tearDown() throws Exception {
		sudoku.clear();
	}

	@Test
	public void testAdd() {
		sudoku.add(1,1,1);
		assertTrue(sudoku.get(1,1) == 1);
		sudoku.add(1,1,2);
		assertTrue(sudoku.get(1,1) == 2);
	}
	
	@Test
	public void testClear() {
		sudoku.add(1,1,1);
		sudoku.add(2,1,2);
		sudoku.add(4,7,4);
		sudoku.add(5,6,9);
		sudoku.clear();
		assertTrue(sudoku.get(1,1) == 0);
		assertTrue(sudoku.get(2,1) == 0);
		assertTrue(sudoku.get(4,7) == 0);
		assertTrue(sudoku.get(5,6) == 0);
	}
	
	@Test
	public void testSolveEmpty() {
		assertTrue(sudoku.solve());
	}

	@Test
	public void testSolveNotEmpty() {
		sudoku.add(0, 2, 8); 
		sudoku.add(0, 5, 9);
		sudoku.add(0, 7, 6);
		sudoku.add(0, 8, 2);
		
		sudoku.add(1, 8, 5);
		
		sudoku.add(2, 0, 1);
		sudoku.add(2, 2, 2);
		sudoku.add(2, 3, 5);
		
		sudoku.add(3, 3, 2);
		sudoku.add(3, 4, 1);
		sudoku.add(3, 7, 9);
		
		sudoku.add(4, 1, 5);
		sudoku.add(4, 6, 6);
		
		sudoku.add(5, 0, 6);
		sudoku.add(5, 7, 2);		
		sudoku.add(5, 8, 8);
		
		sudoku.add(6, 0, 4);
		sudoku.add(6, 1, 1);
		sudoku.add(6, 3, 6);
		sudoku.add(6, 5, 8);
		
		sudoku.add(7, 0, 8);
		sudoku.add(7, 1, 6);
		sudoku.add(7, 4, 3);
		sudoku.add(7, 6, 1);
		
		sudoku.add(8, 6, 4);
		
		assertTrue(sudoku.solve());
	}
	
	@Test
	public void testSolveWrongInput() {
		sudoku.add(0, 2, 8); 
		sudoku.add(0, 5, 9);
		sudoku.add(0, 7, 6);
		sudoku.add(0, 8, 2);
		
		sudoku.add(1, 8, 5);
		
		sudoku.add(2, 0, 1); //här
		sudoku.add(2, 2, 2);
		sudoku.add(2, 3, 5);
		
		sudoku.add(3, 3, 2);
		sudoku.add(3, 4, 1);
		sudoku.add(3, 7, 9);
		
		sudoku.add(4, 1, 5);
		sudoku.add(4, 6, 6);
		
		sudoku.add(5, 0, 6);
		sudoku.add(5, 7, 2);		
		sudoku.add(5, 8, 8);
		
		sudoku.add(6, 0, 4);
		sudoku.add(6, 1, 1);
		sudoku.add(6, 3, 6);
		sudoku.add(6, 5, 8);
		
		sudoku.add(7, 0, 1); //här
		sudoku.add(7, 1, 6);
		sudoku.add(7, 4, 3);
		sudoku.add(7, 6, 1);
		
		sudoku.add(8, 6, 4);
		
		assertFalse(sudoku.solve());
	}
	@Test
	public void testSolveWorldsHardestSudoku() {
		sudoku.add(0, 0, 8);
		
		sudoku.add(1, 2, 3);
		sudoku.add(1, 3, 6);
		
		sudoku.add(2, 1, 7);
		sudoku.add(2, 4, 9);
		sudoku.add(2, 6, 2);
		
		sudoku.add(3, 1, 5);
		sudoku.add(3, 5, 7);
		
		sudoku.add(4, 4, 4);
		sudoku.add(4, 5, 5);
		sudoku.add(4, 6, 7);
		
		sudoku.add(5, 3, 1);
		sudoku.add(5, 7, 3);
		
		sudoku.add(6, 2, 1);
		sudoku.add(6, 7, 6);
		sudoku.add(6, 8, 8);
		
		sudoku.add(7, 2, 8);
		sudoku.add(7, 3, 5);
		sudoku.add(7, 7, 1);
		
		sudoku.add(8, 1, 9);
		sudoku.add(8, 6, 4);
		
		assertTrue(sudoku.solve());		
	}
	
	@Test
	public void testNotSolvableSudoku() {
		
		
		sudoku.add(0, 1, 7);
		sudoku.add(0, 5, 6);
		
		sudoku.add(1, 0, 9);
		sudoku.add(1, 7, 4);
		sudoku.add(1, 8, 1);
		
		sudoku.add(2, 2, 8);
		sudoku.add(2, 5, 9);
		sudoku.add(2, 7, 5);
		
		sudoku.add(3, 1, 9);
		sudoku.add(3, 5, 7);
		sudoku.add(3, 8, 2);
		
		sudoku.add(4, 2, 3);
		sudoku.add(4, 6, 8);
		
		sudoku.add(5, 0, 4);
		sudoku.add(5, 3, 8);
		sudoku.add(5, 7, 1);
		
		sudoku.add(6, 1, 8);
		sudoku.add(6, 3, 3);
		sudoku.add(6, 6, 9);
		
		sudoku.add(7, 0, 1);
		sudoku.add(7, 1, 6);
		sudoku.add(7, 8, 7);
		
		sudoku.add(8, 3, 5);
		sudoku.add(8, 7, 8);
		
		assertFalse(sudoku.solve());
	}
	}
	

