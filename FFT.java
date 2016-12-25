import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.Console;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;



public class FFT {
	
	// convert RGB int matrix into image
	public static BufferedImage RGBtoImage (int[][][]m){ 
		int w = m.length;
		int h = m[0].length;
		BufferedImage output_img = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
		for (int j = 0; j < h; j++){
			for (int i = 0; i < w; i++){
				Color color = new Color(m[i][j][0],m[i][j][1],m[i][j][2]);
				output_img.setRGB(i,j, color.getRGB());
			}
		}
		return output_img;
	}
	
	// One dimensional Fast Fourier Transform, online open source
    public static Complex[] fft1(Complex[] array) {
        int n = array.length;

        // base case
        if (n == 1) return new Complex[] { array[0] };

        // radix 2 Cooley-Tukey FFT

        if (n % 2 != 0) { throw new RuntimeException("n is not a power of 2"); }

        // fft of even terms
        Complex[] even = new Complex[n/2];
        for (int k = 0; k < n/2; k++) {
            even[k] = array[2*k];
        }
        Complex[] q = fft1(even);

        // fft of odd terms
        Complex[] odd  = even;  // reuse the array
        for (int k = 0; k < n/2; k++) {
            odd[k] = array[2*k + 1];
        }
        Complex[] r = fft1(odd);

        // combine
        Complex[] y = new Complex[n];
        for (int k = 0; k < n/2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + n/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }
    
    //Two dimensional FFT
    public static Complex[][][] fft2(Complex[][][] m){
    	int w = m.length;
    	int h = m[0].length;
    	
    	Complex[][][] output = new Complex[w][h][3];
    	Complex [] row = new Complex [w];
    	Complex [] col = new Complex [h];
    	for (int k = 0; k <3; k++){
    		for (int j = 0; j < h; j++){
    			for (int i = 0; i < w; i++){
    				row[i] = m[i][j][k];			
    			}
    			row = fft1(row);
    			for (int i = 0; i < w; i++)
    				output[i][j][k] = row[i];
    		}
    		for (int i = 0; i < w; i++){
    			for (int j = 0; j < h; j++){
    				col[j] = output[i][j][k];			
    			}
    			col = fft1(col);
    			for (int j = 0; j < h; j++)
    				output[i][j][k] = col[j];
    		}
		} 	
    	return output;  	
    }


    // inverse FFT, ifft online open source
    public static Complex[] ifft1(Complex[] x) { 
        int n = x.length;
        Complex[] y = new Complex[n];

        // take conjugate
        for (int i = 0; i < n; i++) {
            y[i] = x[i].conjugate();
        }
        // compute forward FFT
        y = fft1(y);
        // take conjugate again
        for (int i = 0; i < n; i++) {
            y[i] = y[i].conjugate();
        }
        // divide by n
        for (int i = 0; i < n; i++) {
            y[i] = y[i].scale(1.0 / n);
        }
        return y;
    }
    
    //two dimensional inverse fast Fourier Transform 
    public static Complex[][][] ifft2 (Complex[][][] m){
    	int w = m.length;
    	int h = m[0].length;
    	Complex[][][] output = new Complex[w][h][3];
    	Complex [] row = new Complex [w];
    	Complex [] col = new Complex [h];
    	for (int k = 0; k <3; k++){
    		for (int j = 0; j < h; j++){
    			for (int i = 0; i < w; i++){
    				row[i] = m[i][j][k];			
    			}
    			row = ifft1(row);
    			for (int i = 0; i < w; i++)
    				output[i][j][k] = row[i];
    		}
    		for (int i = 0; i < w; i++){
    			for (int j = 0; j < h; j++){
    				col[j] = output[i][j][k];			
    			}
    			col = ifft1(col);
    			for (int j = 0; j < h; j++)
    				output[i][j][k] = col[j];
    		}
		} 	
    	return output;  	
    }

    
    // convert image to matrix, also padding 0 if size is not power of 2 integer
	public static int[][][] imgToMatrix(BufferedImage img, int w, int h){ 
		int w1;
		int h1;
		int[][][] output;
		if ((w & (w - 1)) != 0  || (h & (h - 1)) != 0 ){ // check if it is w or h is power of 2 int
			if ((w & (w - 1)) != 0 ) {w1 = (int)Math.pow(2,Math.ceil(Math.log(w)/Math.log(2)));}
			else{w1 = w;}
			if ((h & (h - 1)) != 0){h1 = (int)Math.pow(2,Math.ceil(Math.log(h)/Math.log(2)));}
			else{h1 = h;}

			output = new int [w1][h1][3];
			for (int k = 0; k <3; k++){
				for (int j = 0; j < h1; j++){
					for (int i = 0; i < w1; i++){
						output[i][j][k] = 0;
					}
				}
			}
		}
		else{
			output = new int [w][h][3];
		}
		
		for (int j = 0; j < h; j++){
			for (int i = 0; i < w; i++){
				int pixel = img.getRGB(i,j);
				Color color = new Color(pixel);
				output[i][j][0] = color.getRed();
				output[i][j][1] = color.getGreen();
				output[i][j][2] = color.getBlue();
			}
		}

		return output;
	}
	
	//The high frequency values will be shifted to the center in spectrum after applying "shift_to_center"
	public static int[][][] shift_to_center(int[][][] matrix){ 
		int w = matrix.length;
		int h = matrix[0].length;
		
		for (int k = 0; k <3; k++){
			for (int j = 0; j < h; j++){
				for (int i = 0; i < w; i++){
					if ((i+j)%2 != 0){
						matrix[i][j][k] = ~matrix[i][j][k];
					}
				}
			}
		}
		return matrix;
	}
	
	//convert from int matrix to complex number matrix 
	public static Complex[][][] getComplexMatrix(int [][][] matrix){ 
		int w = matrix.length;
		int h = matrix[0].length;
		Complex [][][] output = new Complex [w][h][3];
		for (int k = 0; k <3; k++){
			for (int j = 0; j < h; j++){
				for (int i = 0; i < w; i++){
					output[i][j][k] = new Complex((double)(matrix[i][j][k]), (double)(0));
				}
			}
		}
		return output;
	}

	public static int max(int a, int b){
		if (a>=b){
			return a;
		}
		else{
			return b;
		}
	}
	public static int min(int a, int b){
		if (a<b){
			return a;
		}
		else{
			return b;
		}
	}
	public static int scale(int a){ //keep image values within 0 to 255
		if (a<0){
			return 0;
		}
		else if (a>255){
			return 255;
		}
		else{
			return a;
		}
	}
	
	
	// convert complex numbers matrix into image, spectrum = True means we want spectrum from this matrix,
	// otherwise it will be converted into regular image
	public static BufferedImage complexToImg (Complex [][][] matrix, Boolean spectrum){
																						
		int w = matrix.length;
		int h = matrix[0].length;
		BufferedImage output = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
		int [][][] m = new int [w][h][3];
		int max = 0;
		int min = 255;
		for (int k = 0; k <3; k++){
			for (int j = 0; j<h; j++){
				for (int i = 0; i<w; i++){
					if (spectrum){
						m[i][j][k] = (int)Math.log(matrix[i][j][k].abs()+1);
					}
					else{
						m[i][j][k] = (int)matrix[i][j][k].abs();
						m[i][j][k] = scale(m[i][j][k]);
					}
					max = max(max, m[i][j][k]);
				}
			}
		}

		if (spectrum){
			int C = 255/max;	//cofficient C used to scale the image
			for (int k = 0; k <3; k++){
				for (int j = 0; j<h; j++){
					for (int i = 0; i<w; i++){
						m[i][j][k] = (int)(m[i][j][k]*C);				
					}
				}
			}
		}
		
		output = RGBtoImage(m);
		return output;
	}
	
	//Circular filter, ro is cutoff frequency, LowPass = False means highpass
	public static Complex [][] circularFilter(int w, int h, int r0, Boolean LowPass){
		Complex [][] filter = new Complex [w][h];
		int centerX = w/2;
		int centerY = h/2;
		
		for (int j = 0; j < h; j++){
			for (int i = 0; i < w; i++){
				int r = (int)Math.hypot(i-centerX,j-centerY);
				if (LowPass){
					if (r <= r0){
						filter[i][j] = new Complex(1,0);						
					}
					else{
						filter[i][j] = new Complex(0,0);
					}
				}
				else{
					if (r > r0){
						filter[i][j] = new Complex(1,0);		
					}
					else{
						filter[i][j] = new Complex(0,0);
					}
				}			
			}
		}	
		return filter;	
	}
	
	
	//Butterworth filter, ro is cutoff frequency, p is the power, LowPass = False means highpass
	public static Complex [][] ButterWorth(int w, int h, int r0, int p, Boolean LowPass){
		Complex [][] H = new Complex [w][h];
		int centerX = w/2;
		int centerY = h/2;
		for (int j = 0; j < h; j++){
			for (int i = 0; i < w; i++){
				double r = Math.hypot(i-centerX,j-centerY);
				double z = 0;
				if (LowPass){
					z = r/(double)r0;
				}
				else{
					z = (double)r0/r;
				}
				double re = 1/(1+Math.pow(z, 2*p));
				H[i][j] = new Complex(re,0);
			}
		}
		return H;		
	}
	

	//frequency domain image F times frequency domain filter H
	public static Complex [][][] FtimesH(Complex[][][] F, Complex[][] H){ 
		int w = F.length;
		int h = F[0].length;
		Complex [][][] result = new Complex [w][h][3];
		for (int k = 0; k <3; k++){
			for (int j = 0; j < h; j++){
				for (int i = 0; i < w; i++){
					result[i][j][k] = F[i][j][k].times(H[i][j]);
				
				}
			}
		}
		return result;
	}
	
	//w h is the original size of image, Crop the padded image to original size
	public static Complex [][][] crop(Complex [][][] F, int w, int h){ 
		Complex [][][] output = new Complex [w][h][3];
		for (int k = 0; k <3; k++){
			for (int j = 0; j < h; j++){
				for (int i = 0; i < w; i++){
					output[i][j][k] = F[i][j][k];
				
				}
			}
		}
		return output;
	}
	
	
	
	
	//*****************************************Main********************************************************//
	public static void main(String args[])throws IOException{
		// command, '-l' :lowpass; '-h' :highpass; '-c' :circular filter
		// '-b' : butterworth; 
		
		
	    BufferedImage image = null;
	    File file = null;
	    int n = 0;
	    Console console = System.console();
	    String filename = console.readLine("Enter photo filename: ");
  
	    //read image
	    try{
	    	// in java path = "src/images/name.PNG"
	    	// app path = name.PNG
//	    	String image_name = args[0];
	    	String path = new File(filename).getAbsolutePath();
//	    	String path = new File("src/images/Gaussian.jpg").getAbsolutePath();
	    	file = new File(path);
	    	image = ImageIO.read(file);
	      
	    }catch(IOException e){
	      System.out.println(e);
	    }
	    
	    
	    JFrame frame3 = new JFrame();	    //show original image to user
	    frame3.getContentPane().setLayout(new FlowLayout());
	    frame3.getContentPane().add(new JLabel(new ImageIcon(image)));
	    frame3.pack();
	    frame3.setVisible(true);
	    
	    
	    
	    //**********************display original spectrum of image, and perform FFT ***********************//
	    int w0 = image.getWidth();
		int h0 = image.getHeight();
		System.out.println("Image size is " + w0 + " x " + h0);
		System.out.println("Total Pixels: " + w0*h0);
		long startTime = System.currentTimeMillis();
	    int [][][] colors = shift_to_center(imgToMatrix(image, w0, h0));
		Complex [][][] f = getComplexMatrix(colors);
		Complex [][][] F = fft2(f);
		
		
		BufferedImage img_spectrum = complexToImg(F, true);
	    JFrame frame = new JFrame();	    //display spectrum of original image to user
	    frame.getContentPane().setLayout(new FlowLayout());
	    frame.getContentPane().add(new JLabel(new ImageIcon(img_spectrum)));
	    frame.pack();
	    frame.setVisible(true);
	    
	    long endTime   = System.currentTimeMillis();
	    long totalTime = endTime - startTime;
	    System.out.println("FFT2 Run Time is: ");
	    System.out.println(totalTime + " ms");
	    
	    int w = F.length;
		int h = F[0].length;
		Complex [][] H = new Complex [w][h]; 
	    
		
		
	    
	    //********************** Now Create filter by user**********************************//
	    System.out.println("Please choose filter type, -l is lowpass, -h is highpass");
	    String type = console.readLine("Enter filter type: ");
	    System.out.println("Now pick a filter method, -c for circular, -b for butterworth");
	    String method = console.readLine("Enter filter method: ");
	    
		System.out.printf("Pick cutoff radius r0, from range 0 to %d", (int)(max(w,h)/2));
		System.out.println();
		String R0 = console.readLine("Specify r0: ");
		int r0 = Integer.parseInt(R0);
		if (method.equals("-b")){
			System.out.println("Choose the order of filter n, from 0, 1, 2...");
			String power = console.readLine("Specify n: ");
			n = Integer.parseInt(power);
	    }	
		
		if (method.equals("-b")){
			if (type.equals("-l")){
				H = ButterWorth(w,h, r0, n, true);
			}
			else if (type.equals("-h")){
				H = ButterWorth(w,h, r0, n, false);
			}
			else{
				System.out.println("Error: wrong filter type input");
				System.exit(0);
			}
		}
		else if (method.equals("-c")){
			if (type.equals("-l")){
				H = circularFilter(w,h, r0, true);
			}
			else if (type.equals("-h")){
				H = circularFilter(w,h, r0, false);
			}
			else{
				System.out.println("Error: wrong filter type input");
				System.exit(0);
			}
		}
		else{
			System.out.println("Error: wrong filter method input");
			System.exit(0);
		}

		Complex [][][] m_filt = FtimesH(F, H); //apply filter on image F

		BufferedImage fftImg = complexToImg(m_filt,true); //filtered spectrum		
		Complex [][][] filted = ifft2(m_filt);
		
		if (w != w0 || h != h0) filted = crop(filted, w0, h0); //if size changed, then resize it
															// w0 and h0 is the original size of image
		BufferedImage fftImge = complexToImg(filted,false); //filtered image
		
		String path = new File(".").getAbsolutePath();
	    path = path.substring(0, path.length()-1);
	    ImageIO.write(fftImg, "png", new File(path + "fftFiltImg.jpg")); //save filtered spectrum image to current folder
	    System.out.println("File fftFiltImg.jpg created -- spectrum of filtered image");
	    ImageIO.write(fftImge, "png", new File(path + "filtImage.jpg")); //save filtered image to current folder
	    System.out.println("File filtImage.jpg created -- filtered image");
	    

	    JFrame frame1 = new JFrame();	    
	    frame1.getContentPane().setLayout(new FlowLayout());
	    frame1.getContentPane().add(new JLabel(new ImageIcon(fftImg)));
	    frame1.pack();
	    frame1.setVisible(true);
	    
	    JFrame frame2 = new JFrame();	    
	    frame2.getContentPane().setLayout(new FlowLayout());
	    frame2.getContentPane().add(new JLabel(new ImageIcon(fftImge)));
	    frame2.pack();
	    frame2.setVisible(true);
	    System.out.println("Press 'ctrl + c' to quit program");
	}
}
