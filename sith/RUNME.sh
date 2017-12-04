if [ -d output/ ]; then
   rm -rf output/
   echo "Removed the directory output and now I am going to create it again."
else
   echo "Directory output doesn't exist. Going to create it"
fi

mkdir output/
cp input/* output
cd output
echo "These are the parameters that are going to be used for the taboo search:"
cat sith_input.txt
read -p  "Please check the doc I haven't written yet for more information. If you are happy with this parameters, please continue. Otherwise, modify whatever you think needs to be modified and run the script again. Do you want to continue? [y/n]
" choice 

while [[ $choice != "y" && $choice != "n" ]]; do
  read -p "I said, do you want to continue? [y/n]
  " choice
done

if [ $choice == "y" ]; then
   python sith.py -i sith_input.txt
else
   echo "ABORTING"
   exit
fi
