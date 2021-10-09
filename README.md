# Construction of Active and Passive circuits

This repository contains the implementation of four filters in Multisim:

- Butterworth Low-Pass
- Chebyshev High-Pass
- Inverse Chebyshev Band-Pass
- Chebyshev Band-Elimination

Moreover in each folder you can find the MATLAB scripts for the calculation of the parameters of the electronic components.



## Butterworth Low-Pass

After the parameters calculation the following transfer function results:



<p allign='center'>
    <img src='images/lp_tf.jpg' width='60%'>
</p>



This function is graphically illustrated in MATLAB



<p allign='center'>
    <img src='images/lp_tf_plot.jpg' width='49%'> <img src='images/lp_de.jpg' width='49%'>
</p>



Calculating the values of the electronic components, the following circuit is designed:

<p allign='center'>
    <img src='images/lp_circuit.jpg'>
</p>



To investigate the proper function of this circuit, a pulse was given as input:

<p allign='center'>
    <img src='images/lp_input.jpg' width='60%'>
</p>

Subsequently an oscillator was used in the input and the output of the circuit:

<p allign='center'>
    <img src='images/lp_circuit_2.jpg'>
</p>

 The transient analysis using this oscillator resulted in the following graph (red: input, green: output):

<p allign='center'>
    <img src='images/lp_input_output.jpg'>
</p>

Finally, a Fourier analysis was performed for the input and the output signal:

<p allign='center'>
    <img src='images/lp_input_fourier.jpg'>
</p>



<p allign='center'>
    <img src='images/lp_output_fourier.jpg'>
</p>



From this graphs we can conclude that the circuit is working properly. The high frequencies are not allowed to pass through the filter.

---





## Chebyshev High-Pass

After the parameters calculation the following transfer function results:



<p allign='center'>
    <img src='images/hp_tf.jpg' width='60%'>
</p>



This function is graphically illustrated in MATLAB



<p allign='center'>
    <img src='images/hp_tf_plot.jpg' width='49%'> <img src='images/hp_de.jpg' width='49%'>
</p>



Calculating the values of the electronic components, the following circuit is designed:

<p allign='center'>
    <img src='images/hp_circuit.jpg'>
</p>



To investigate the proper function of this circuit, a sum of sine signals was given as input and an oscillator was used:

<p allign='center'>
    <img src='images/hp_circuit_2.jpg'>
</p>

 The transient analysis using this oscillator resulted in the following graph for input:

<p allign='center'>
    <img src='images/hp_input_signal.jpg'>
</p>

and the following graph for the output signal:

<p allign='center'>
    <img src='images/hp_output_signal.jpg'>
</p>



Finally, a Fourier analysis was performed for the input and the output signal:

<p allign='center'>
    <img src='images/hp_input_fourier.jpg'>
</p>



<p allign='center'>
    <img src='images/hp_output_fourier.jpg'>
</p>



From this graphs we can conclude that the circuit is working properly. The low frequencies are not allowed to pass through the filter.

---



## Inverse Chebyshev Band-Pass

After the parameters calculation the following transfer function results:

<p allign='center'>
    <img src='images/bp_tf.jpg' width='80%'>
</p>



This function is graphically illustrated in MATLAB



<p allign='center'>
    <img src='images/bp_tf_plot.jpg' width='49%'> <img src='images/bp_de.jpg' width='49%'>
</p>



Calculating the values of the electronic components, the following circuit is designed:

<p allign='center'>
    <img src='images/bp_circuit.jpg'>
</p>



To investigate the proper function of this circuit, a sum of sine signals was given as input and an oscillator was used:

<p allign='center'>
    <img src='images/bp_circuit_2.jpg'>
</p>

 The transient analysis using this oscillator resulted in the following graph for input (red) and output (green):

<p allign='center'>
    <img src='images/bp_input_output.jpg'>
</p>



Finally, a Fourier analysis was performed for the input and the output:

<p allign='center'>
    <img src='images/bp_input_fourier.jpg'>
</p>



<p allign='center'>
    <img src='images/bp_output_fourier.jpg'>
</p>



From this graphs we can conclude that the circuit is working properly. The low and the high frequencies are not allowed to pass through the filter.

---



## Chebyshev Band-Elimination

After the parameters calculation the following transfer function results:

<p allign='center'>
    <img src='images/be_tf.jpg' width='80%'>
</p>



This function is graphically illustrated in MATLAB



<p allign='center'>
    <img src='images/be_tf_plot.jpg' width='49%'> <img src='images/be_de.jpg' width='49%'>
</p>



Calculating the values of the electronic components, the following circuit is designed:

<p allign='center'>
    <img src='images/be_circuit.jpg'>
</p>



To investigate the proper function of this circuit, a sum of sine signals was given as input and an oscillator was used:

<p allign='center'>
    <img src='images/be_circuit_2.jpg'>
</p>

 The transient analysis using this oscillator resulted in the following graph for input and output:

<p allign='center'>
    <img src='images/be_input_signal.jpg'>
</p>



<p allign='center'>
    <img src='images/be_output_signal.jpg'>
</p>



Finally, a Fourier analysis was performed for the input and the output:

<p allign='center'>
    <img src='images/be_input_fourier.jpg'>
</p>



<p allign='center'>
    <img src='images/be_output_fourier.jpg'>
</p>



From this graphs we can conclude that the circuit is working properly. The in-between frequencies (in this case 2.287 𝑘𝐻𝑧) are not allowed to pass through the filter.

---

