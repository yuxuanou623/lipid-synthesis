import requests

# Replace this with your actual Flask server URL and port
url = 'http://localhost:8080/reaction'

# Example payload: a list of dictionaries
data = [
    {
        "head": "CC#CCCCCCCCCCN",
        "tail_smiles1": "COCCN(CC(=O)O)C[C@H]1C[C@@H](CO)CO1"
    },
    {
        "head": "CC#CCCCCCCCCCN",
        "tail_smiles1": "CN(C)CCN1CCN(C(=O)[C@H](N)CO)CC1"
    },
    # Add more dictionaries as needed
]

response = requests.post(url, json=data)

# Print the server's response
print("Status code:", response.status_code)
print("Response JSON:", response.json())
