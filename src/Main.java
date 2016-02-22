/*-----------------------------------------------------------------
 * CS191 Major Output 2
 * Written by: Mestica Samela G. Papna
 * -----------------------------------------------------------------
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import java.util.LinkedList;
import java.util.Collections;
import java.util.Random;
import java.util.Arrays;

import javax.swing.JFrame;
import java.awt.EventQueue;
import org.math.plot.*;

import javax.swing.JPanel;
import java.awt.Color;



/* CHROMOSOME OBJECT
 * 	Basic Object for this application that uses Evolutionary Algorithm.
 */
class Chromosome {
	double beta, gamma;
   RMSD rmsd;
	
	public void setParameters(double beta, double gamma) {
		this.beta = beta;
		this.gamma = gamma;
	} 
	
	public double beta() {
		return this.beta;
	} 

	public double gamma() {
		return this.gamma;
	} 
	
	public void setRMSD(RMSD rmsd) {
		this.rmsd = rmsd;
	} 
	
	public RMSD getRMSD() {
		return this.rmsd;
	} 
	
	public double getFitnessValue() {
		return (this.rmsd.getRMSDS() + this.rmsd.getRMSDI() + this.rmsd.getRMSDR()) / 3.0;
	}
}

/* SIR OBJECT
 * 	Creation of the data type that will hold the set of values of S_generated, 
 * 	I_generated and R_generated plus the S_actual, I_actual and R_actual that
 * 	are used in the computation found in the EpidemicModelling Class. 
 */
class SIR {
	Double[] S;
	Double[] I;
	Double[] R;
	
	public void setParameters(Double[] S,Double[] I, Double[] R) {
		this.S = S;
		this.I = I;
		this.R = R;		 
	}
	
	public Double[] getS() {
		return this.S;
	}
	
	public Double[] getI() {
		return this.I;
	}
	
	public Double[] getR() {
		return this.R;
	}
}

/*	RMSD OBJECT
 * 	To be used in the chromosome object, like an enum in C, to store the values of the RMSD(S(t)), RMSD(I(t))
 * 	and RMSD(R(t)) of the chromosomes in each generation
 */
class RMSD {
	double rmsdS, rmsdI, rmsdR;

   public void setParameters(double rmsdS, double rmsdI, double rmsdR) {
   	this.rmsdS = rmsdS;
   	this.rmsdI = rmsdI;
   	this.rmsdR = rmsdR;
   }
   
   public double getRMSDS() {
   	return this.rmsdS;
   }

   public double getRMSDI() {
   	return this.rmsdI;
   }

   public double getRMSDR() {
   	return this.rmsdR;
   }
}

/* EPIDEMIC MODELLING APPLICATION
 * 	Glues the three objects defined above to implement 
 *    the application specified for CS191 Major Output 2
 */
class EpidemicModelling {
	double F = 1.0; // Static mutation factor F
	double CR = 0.8; // 80% Crossover Ratio, this is at somepoint useless since I am using exp for crossover but nevertheless I was alotted with it 
	
	
	public double k1S(double b, double n, double s, double i) {
		return -(b/n)*s*i;
	}

	public double k1I(double b, double n, double s, double i, double g) {
		return (b/n)*s*i - g*i;
	}
	
	public double k1R(double g, double i) {
		return g*i;
	}
	
	public double k2S(double b, double n, double s, double i, double h, double k1S, double k1I) {
		return -(b/n)*(s+h*0.5*k1S)*(i+h*0.5*k1I);
	}
	
	public double k2I(double b, double n, double s, double i, double g, double h, double k1S, double k1I) {
		return (b/n)*(s+h*0.5*k1S)*(i+h*0.5*k1I) - g*(i+h*0.5*k1I);
	}

	public double k2R(double g, double i, double k1I, double h) {
		return g*(i+h*0.5*k1I);
	}
	
	public double k3S(double b, double n, double s, double i, double h, double k2S, double k2I) {
		return -(b/n)*(s+h*0.5*k2S)*(i+h*0.5*k2I);
	}
	
	public double k3I(double b, double n, double s, double i, double g, double h, double k2S, double k2I) {
		return (b/n)*(s+h*0.5*k2S)*(i+h*0.5*k2I) - g*(i+h*0.5*k2I);
	}

	public double k3R(double g, double i, double k2I, double h) {
		return g*(i+h*0.5*k2I);
	}
	
	public double k4S(double b, double n, double s, double i, double h, double k3S, double k3I) {
		return -(b/n)*(s+h*k3S)*(i+h*k3I);
	}
	
	public double k4I(double b, double n, double s, double i, double g, double h, double k3S, double k3I) {
		return (b/n)*(s+h*k3S)*(i+h*k3I) - g*(i+h*k3I);
	}

	public double k4R(double g, double i, double k3I, double h) {
		return g*(i+h*k3I);
	}


	public SIR rk4(double beta, double gamma, double initialS, double initialI, double initialR, int steps) {
		
		double N = initialI + initialR + initialS;
		
		Double[] S = new Double[steps];
		Double[] I = new Double[steps];
		Double[] R = new Double[steps];
	
		double k1S, k2S, k3S, k4S;
		double k1I, k2I, k3I, k4I;
		double k1R, k2R, k3R, k4R;
				
		/*
		 * We will find the value for the first set of SIR 
		 * that is S_1, I_1, R_1 but termed as S[0], I[0], R
		 * [0], in here. 
		 */		
		k1S = k1S(beta, N, initialS, initialI);		
		k1I = k1I(beta, N, initialS, initialI, gamma);
		k1R = k1R(gamma, initialI);
		
		k2S = k2S(beta, N, initialS, initialI, 1, k1S, k1I);
		k2I = k2I(beta, N, initialS, initialI, gamma, 1, k1S, k1I);
		k2R = k2R(gamma, initialI, k1I, 1);
		
		k3S = k3S(beta, N, initialS, initialI, 1, k2S, k2I);
		k3I = k3I(beta, N, initialS, initialI, gamma, 1, k2S, k2I);
		k3R = k3R(gamma, initialI, k2I, 1);

		k4S = k4S(beta, N, initialS, initialI, 1, k3S, k3I);
		k4I = k4I(beta, N, initialS, initialI, gamma, 1, k3S, k3I);
		k4R = k4R(gamma, initialI, k3I, 1);
		
		/*
		 * After computing all the Ks that we need.
		 */
		S[0] = initialS + (1/6.0)*(k1S + 2*k2S + 2*k3S + k4S);		
		I[0] = initialI + (1/6.0)*(k1I + 2*k2I + 2*k3I + k4I);		
		R[0] = initialR + (1/6.0)*(k1R + 2*k2R + 2*k3R + k4R);
		
				
		/*
		 * Compute for the rest of the sets of S, I and R.
		 * We will also do what we have done on the initial
		 * values of SIR
		 */
		for(int k = 1; k < steps; k++) {
			double h = 1.0; //Defining h as the step size
			
			k1S = k1S(beta, N, S[k-1], I[k-1]);		
			k1I = k1I(beta, N, S[k-1], I[k-1], gamma);
			k1R = k1R(gamma, I[k-1]);
		
			k2S = k2S(beta, N, S[k-1], I[k-1], h, k1S, k1I);
			k2I = k2I(beta, N, S[k-1], I[k-1], gamma, h, k1S, k1I);
			k2R = k2R(gamma, I[k-1], k1I, h);
		
			k3S = k3S(beta, N, S[k-1], I[k-1], h, k2S, k2I);
			k3I = k3I(beta, N, S[k-1], I[k-1], gamma, h, k2S, k2I);
			k3R = k3R(gamma, I[k-1], k2I, h);

			k4S = k4S(beta, N, S[k-1], I[k-1], h, k3S, k3I);
			k4I = k4I(beta, N, S[k-1], I[k-1], gamma, h, k3S, k3I);
			k4R = k4R(gamma, I[k-1], k3I, h);

			S[k] = S[k-1] + (h/6.0)*(k1S + 2*k2S + 2*k3S + k4S);		
			I[k] = I[k-1] + (h/6.0)*(k1I + 2*k2I + 2*k3I + k4I);		
			R[k] = R[k-1] + (h/6.0)*(k1R + 2*k2R + 2*k3R + k4R);			
		}  
			
		SIR temp = new SIR();
		temp.setParameters(S, I, R);
	
		return temp;
	}
	
	/*
	 * Gathers the values from the SIRStats file
	 */
	public SIR readData() {
		SIR temp = new SIR();
		
		try {
			FileReader inSIRStats = new FileReader("SIRStats");
			BufferedReader bufferReaderInSIRStats = new BufferedReader(inSIRStats);
			
			List<String> elephantList;
			ArrayList<ArrayList <Double>> penguinList = new ArrayList<ArrayList <Double>>();
			
			int steps = 0;
			String lineCurrent;

			while((lineCurrent = bufferReaderInSIRStats.readLine()) != null) {
				penguinList.add(new ArrayList<Double>());

				elephantList = Arrays.asList(lineCurrent.split("\\s+"));				
				for(int i = 0; i < elephantList.size(); i++) {
					penguinList.get(steps).add(Double.parseDouble(elephantList.get(i)));
				}
				steps++;				
			} 

			Double[] S = new Double[steps];
			Double[] I = new Double[steps];
			Double[] R = new Double[steps];			

			for(int k = 0; k < steps; k++) {
				S[k] = penguinList.get(k).get(0);
				I[k] = penguinList.get(k).get(1);
				R[k] = penguinList.get(k).get(2);				 
			}
			
			temp.setParameters(S, I, R);
		} catch(Exception e) {
			System.out.println("readData(): Error while reading file line by line");			
		}
		
		return temp;
	}
	
	/*
	 * The function that computes the rmsd for each chromosome and its corresponding 
	 * beta and gamma. 
	 */
	public RMSD rmsd(double beta, double gamma, SIR generated, SIR actual) {
		RMSD temp = new RMSD();

		try {
			double N = actual.getS()[0] + actual.getI()[0] + actual.getR()[0];
			double summation = 0;
	
			/*
			 * We have set the length of the S[], I[], R[] in generated
			 * as the bound in the loop because it has a smaller points than the 
			 * actual, 45 versus 46 respectively, since we have ommitted the S[0], I[0],
			 * R[0] from the given data when using the RK4.
			 */
			 
	
			for(int k = 0; k < generated.getS().length; k++) {
				summation += actual.getS()[k+1] - generated.getS()[k];    // k+1 for actual S because we don't read the initial values of SIR
 			}

			double S = Math.sqrt(summation*summation/N);
			
			summation = 0;

			for(int k = 0; k < generated.getI().length; k++) {
				summation += actual.getI()[k+1] - generated.getI()[k];		// k+1 for actual I because we don't read the initial values of SIR
			}
		
			double I = Math.sqrt(summation*summation/N);
		
			summation = 0;
		
			for(int k = 0; k < generated.getR().length; k++) {
				summation += actual.getR()[k+1] - generated.getR()[k];		// k+1 for actual R because we don't read the initial values of SIR
			}
		
			double R = Math.sqrt(summation*summation/N);

			temp.setParameters(S, I, R);
		} catch(Exception e) {
			System.out.println("rmsd(): " + e);
		}
		
		
		return temp;
	}
		
	/*
	 * Chooses the chromosome with the lowest fitnessValue or average mean error
	 */	
	public Chromosome findBestChromosome(ArrayList<Chromosome> chromosomes) {
		Chromosome bestChromosome = chromosomes.get(0);
		
		for(int k = 1; k < 40; k++) {
			if(bestChromosome.getFitnessValue() > chromosomes.get(k).getFitnessValue()) {
				bestChromosome = chromosomes.get(k);
			}
		}
		
		return bestChromosome;
	}
	
	/*
	 * Start of the Differential Evolution process applied to SIR Epidemiology Model
	 */
	public void startProcess() {
		Random rand = new Random();
		
		try {
			ArrayList<Chromosome> chromosomes = new ArrayList<Chromosome>();

			//INITIALIZATION
			for(int k = 0; k < 40; k++) {
				double beta = rand.nextDouble();
				double gamma = rand.nextDouble();
				chromosomes.add(new Chromosome());
				chromosomes.get(k).setParameters(beta, gamma);
			}
				
			//EVALUATION OF PARENTS
			SIR actual = readData();
			for(int k = 0; k < 40; k++) {
				SIR generated = rk4(chromosomes.get(k).beta(), chromosomes.get(k).gamma(), actual.getS()[0], actual.getI()[0], actual.getR()[0], actual.getR().length - 1);								
				chromosomes.get(k).setRMSD(rmsd(chromosomes.get(k).beta(), chromosomes.get(k).gamma(), generated, actual));
			}
			
			// DE/best/2/exp
			int n; //number of iterations
			int n2 = 0; //number of times the best chromosome has been repeated
			Chromosome[] arrayOfBest = new Chromosome[1000000]; //in each of the indices of this array contains the bestChromosome for its generation
			for(n = 0; n < 1000000 && n2 < 1000; n++) {
				ArrayList<Chromosome> mutants = new ArrayList<Chromosome>();	
				
				//Find the chromosome with the best fitness value (i.e. best = lowest).			
				Chromosome temp = findBestChromosome(chromosomes);

				/* GENERATE MUTANT POPULATION
				 * 	We were given a best/2 approach. As seen from the line of code above we find the best chromosome first
				 *    then randomized for chromosome from the current generation's parent set and apply the formula
				 *    for each parameter(i.e. beta and gamma) which is 
				 */
				for(int k = 0; k < 40; k++) { 	
					Chromosome temp2 = chromosomes.get(rand.nextInt(40));
					Chromosome temp3 = chromosomes.get(rand.nextInt(40));
					Chromosome temp4 = chromosomes.get(rand.nextInt(40));
					Chromosome temp5 = chromosomes.get(rand.nextInt(40));
					
					mutants.add(new Chromosome());
					mutants.get(k).setParameters(temp.beta() + F*(temp2.beta() - temp3.beta() + temp4.beta() - temp5.beta()), temp.gamma() + F*(temp2.gamma() - temp3.gamma() + temp4.gamma() - temp5.gamma())); 

					if(mutants.get(k).beta() <= 0) {
						mutants.get(k).setParameters(0.01, mutants.get(k).gamma());
					}
					if(mutants.get(k).beta() >= 1) {
						mutants.get(k).setParameters(0.99, mutants.get(k).gamma());
					}

					if(mutants.get(k).gamma() <= 0) {
						mutants.get(k).setParameters(mutants.get(k).beta(), 0.01);
					}
					if(mutants.get(k).gamma() >= 1) {
						mutants.get(k).setParameters(mutants.get(k).beta(), 0.99);
					}
				}
				
				/* CROSSOVER
				 * 	Since we were given an exp approach. We randomized two values for beta and gamma.
				 * 	Check for their parity(i.e. even or odd) and with this the beta and gamma values 
				 * 	of the mutant population obtained are subject to change. They either take the beta 
				 * 	or gamma of their parent or retain it.  
				 */ 
				for(int k = 0; k < 40; k++) {
					int randomNumberForBeta = rand.nextInt(100);
					int randomNumberForGamma = rand.nextInt(100);

					if((randomNumberForBeta % 2 == 1) && (randomNumberForGamma % 2 == 1)) {
						mutants.get(k).setParameters(chromosomes.get(k).beta(), mutants.get(k).gamma());
					} else if ((randomNumberForBeta % 2 == 0) && (randomNumberForGamma % 2 == 0)) {
						mutants.get(k).setParameters(mutants.get(k).beta(), chromosomes.get(k).gamma());
					} else if ((randomNumberForBeta % 2 == 1) && (randomNumberForGamma % 2 == 0))	{
						mutants.get(k).setParameters(chromosomes.get(k).beta(), chromosomes.get(k).gamma());
					}		
				}
				
				/* EVALUATION OF MUTANTS
				 * 	Evaluates the mutants to set their fitness values
				 */
				for(int k = 0; k < 40; k++) {
					SIR generated = rk4(mutants.get(k).beta(), mutants.get(k).gamma(), actual.getS()[0], actual.getI()[0], actual.getR()[0], actual.getR().length - 1);								
					mutants.get(k).setRMSD(rmsd(mutants.get(k).beta(), mutants.get(k).gamma(), generated, actual));
				}
				
				/* SELECTION
				 * 	Selects between the mutant and the parents on which among them has the lesser average mean error/ fitness value.
				 */
				for(int k = 0; k < 40; k++) {
					if(mutants.get(k).getFitnessValue() < chromosomes.get(k).getFitnessValue()) {
						chromosomes.get(k).setParameters(mutants.get(k).beta(), mutants.get(k).gamma());
					}
				}

				/* EVALUATION AFTER SELECTION
				 *   A must because the RMSD of the chromosomes will be set to null after 
				 *   the selection process. The chromosomes', which is the new parent here,
				 *   parameters have been reset.
				 */
				for(int k = 0; k < 40; k++) {
					SIR generated = rk4(chromosomes.get(k).beta(), chromosomes.get(k).gamma(), actual.getS()[0], actual.getI()[0], actual.getR()[0], actual.getR().length -1);								
					chromosomes.get(k).setRMSD(rmsd(chromosomes.get(k).beta(), chromosomes.get(k).gamma(), generated, actual));
				}				
								
				/* FIND THE BEST CHROMOSOME & TERMINATION CONDITION TWO UPDATER
				 * 	Finds the best chromosome for the current generation and the 
				 * 	number of times it remained to be so per iteration in the loop.
				 *    It updates the value of n2 which is a variable for terminating
				 * 	condition number two.
				 */
				arrayOfBest[n] = findBestChromosome(chromosomes);				
				if(n > 0) {
					if(arrayOfBest[n].equals(arrayOfBest[n-1])) {
						n2++;
					} else {
						n2 = 0;
					}
				}
				 
				/* TERMINATION CONDITION THREE
				 * 	Ends the iteration once a best chromosome reaches an aveage RMSD less than 0.000001
				 */ 			
				if(arrayOfBest[n].getFitnessValue() < 0.000001) {
					break;
				}
			}
			

			/* FINAL REPORT
			 * 	Prints in the console the specified output for this major output.
			 */
			System.out.println("Number of generations: " + n +"\n\n");
			Chromosome bestChromosome = findBestChromosome(chromosomes);
			System.out.println("The Systems of Equations Using the Best Beta and Gamma");
			double N = actual.getS()[0] + actual.getI()[0] + actual.getR()[0];
			System.out.println("dS/dt = -(" + bestChromosome.beta() +"/"+N+")SI");
			System.out.println("dI/dt = (" + bestChromosome.beta() +"/"+N+")SI - "+bestChromosome.gamma()+"I");
			System.out.println("dR/dt = "+ bestChromosome.gamma() +"I\n\n");
			System.out.println("Final Set of Chromosomes Together with their Beta, Gamma, RMSD(S(t)),  RMSD(I(t)) and RMSD(R(t))");
			for(int k = 0; k < 40; k++) {
				System.out.println("Chromosome[" + (k + 1) + "]: beta= "+chromosomes.get(k).beta()+ " gamma= "+chromosomes.get(k).gamma()+" RMSD(S(t))= " +chromosomes.get(k).getRMSD().getRMSDS()+" RMSD(I(t))= "+chromosomes.get(k).getRMSD().getRMSDI()+" RMSD(R(t))= "+chromosomes.get(k).getRMSD().getRMSDR());
			}
			
			/* PLOT THE <SIR> OF BEST CHROMOSOME
			 * 	Implements the graph for the SIR of the best chromosome found using DE Algorithm
			 * 	and some libraries from jmathplot.
			 */
			JFrame frame = new JFrame("THE SIR MODEL");
			frame.setResizable(true);
			frame.setBounds(100, 100, 891, 595);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			frame.getContentPane().setLayout(null);

			SIR generated = rk4(bestChromosome.beta(), bestChromosome.gamma(), actual.getS()[0], actual.getI()[0], actual.getR()[0], actual.getR().length - 1);								
			double[]	 t = new double[generated.getS().length];
			double[]	 s = new double[generated.getS().length];
			double[]	 i = new double[generated.getS().length];
			double[]	 r = new double[generated.getS().length];
			
			for(int k = 0; k < generated.getS().length;k++) {
				t[k] = k;
				s[k] = (double) generated.getS()[k];
				i[k] = (double) generated.getI()[k];
				r[k] = (double) generated.getR()[k];
			}
			
			JPanel panel = new Plot2DPanel();
			panel.setBounds(15, 11, 870, 515);
			frame.getContentPane().add(panel);

			((Plot2DPanel) panel).addLinePlot("S", Color.blue, t,s);
			((Plot2DPanel) panel).addLinePlot("I", Color.green, t, i);
			((Plot2DPanel) panel).addLinePlot("R", Color.red, t, r);

			frame.setVisible(true);
	
	} catch(Exception e) {
			System.out.println("startProcess(): " + e);
		} finally {
			
		}						
	}
	
	public static void main(String[] args) throws IOException {
		EpidemicModelling em = new EpidemicModelling();
		em.startProcess();	

		/* SHOW THE GRAPH
		 * 	Necessary to in running the SIR Model Graph
		 */
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					EpidemicModelling window = new EpidemicModelling();
					
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
		

	}
}


public class Main {
	public static void main(String args[]) throws IOException {		
		EpidemicModelling.main(args);
	}
} 


