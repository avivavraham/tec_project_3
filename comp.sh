#sudo chmod 755 'filename'
#ssh -L 8080:yonatanromm@gate.tau.ac.il:22 yonatanromm@nova.cs.tau.ac.il   

#python3 setup.py build_ext --inplace
#./spkmeans 3 lnorm input_1.txt
#gcc -ansi -Wall -Wextra -Werror -pedantic-errors main.c -lm -o spkmeans
#python3 kmeans_pp.py 3 100 0.01 test_data/input_1_db_1.txt test_data/input_1_db_2.txt
#python3 kmeans_pp.py 7 0.01 test_data/input_2_db_1.txt test_data/input_2_db_2.txt
#python3 kmeans_pp.py 15 0.01 test_data/input_3_db_1.txt test_data/input_3_db_2.txt
#
#rm -r build
#rm mykmeanssp.cpython-39-darwin.so

gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c nsclustering.c helper.c -lm -o spkmeans
./spkmeans wam /Users/I555060/University/software-project/tec_project_3/testfiles/spk_0.txt &> output.txt
#cat output.txt
#./spkmeans wam testfiles/spk_9.txt | diff - testfiles/outputs/c/wam/spk_9.txt
#./yuval_tester.sh /Users/I555060/University/software-project/tec_project_3/testfiles