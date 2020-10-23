input="IMOS_-_Argo_Profiles_URLs.txt"
while IFS= read -r line
do
  wget -v $line
done < "$input"
