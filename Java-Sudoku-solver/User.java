package Sudoku;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.GridLayout;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

public class User {
	private Sudoku sudoku;
	private Map<Integer, JTextField> map;

	public User(Sudoku sudoku) {
		this.sudoku = sudoku;
		this.map = new HashMap<>();
		SwingUtilities.invokeLater(() -> createWindow(sudoku, "Sudoku", 100, 300));
	}

	private void createWindow(Sudoku sudoku, String title, int width, int height) {
		JFrame frame = new JFrame(title);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container pane = frame.getContentPane();

		// skapar sudokuplanen och knapparna
		JPanel panelup = new JPanel();
		panelup.setLayout(new GridLayout(9, 9));
		for (int i = 0; i < 81; i++) {
			JTextField temp = new JTextField(null);
			int a = i / 3;
			if (a == 1 || a == 4 || a == 7 || a == 9 || a == 11 || a == 12 || a == 14 || a == 15 || a == 17 || a == 19
					|| a == 22 || a == 25) {
				temp.setBackground(Color.gray);
			}
			panelup.add(temp);
			map.put(i, temp);
		}

		// knapparna
		JPanel paneldown = new JPanel();
		JButton solve = new JButton("Solve");
		JButton clear = new JButton("Clear");
		solve.addActionListener(event -> {
			boolean b1 = addNumbersToMatrix();
			if (b1 == false) {
				JOptionPane.showMessageDialog(frame, "Felaktigt tecken");
			} else {
				boolean b2 = sudoku.solve();
				if (b2) {
					showSolve();
				} else {
					JOptionPane.showMessageDialog(frame, "Sudokut har ingen lösning");
					sudoku.clear();
				}
			}
		});

		clear.addActionListener(event -> {
			sudoku.clear();
			showSolve();
		});
		paneldown.add(solve);
		paneldown.add(clear);

		// sätter ihop allt
		pane.add(paneldown, BorderLayout.SOUTH);
		pane.add(panelup, BorderLayout.CENTER);

		frame.pack();
		frame.setVisible(true);

	}

	public boolean addNumbersToMatrix() {
		int r = 0;
		int c = -1;
		for (int i = 0; i < 81; i++) {
			String string = map.get(i).getText();
			if (string.equals("")) {
				c = c + 1;
				if (c == 9) {
					c = 0;
					r++;
				}
			} else {
				try {
					int temp = Integer.valueOf(string.trim());
					if (temp > 0 && temp < 10) { // ändrar från -1
						c = c + 1;
						if (c == 9) {
							c = 0;
							r++;
						}
						sudoku.add(r, c, temp);
					} else {
						return false;
					}
				} catch (NumberFormatException e) {
					return false;
				}
			}
		}
		return true;
	}

	public void showSolve() {
		int r = 0;
		int c = 0;
		for (int i = 0; i < 81; i++) {
			String nbr = String.valueOf(sudoku.get(r, c));
			if (sudoku.get(r, c) == 0) {
				map.get(i).setText("");
			} else {
				map.get(i).setText(nbr);
			}
			c = c + 1;
			if (c == 9) {
				c = 0;
				r++;
			}
		}
	}

}