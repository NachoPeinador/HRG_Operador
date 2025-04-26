# HRG_Operador.sage  
# Versión 2.0 con conexión a LMFDB API  
# Licencia MIT: https://opensource.org/licenses/MIT  

from sage.all import *  
import requests  
import pandas as pd  
import json  

def obtener_ceros_LMFDB(q, chi_label, max_zeros=10):  
    r"""  
    Obtiene ceros no triviales de L(s,χ) desde LMFDB.  

    Parámetros:  
        q (int): Módulo  
        chi_label (str): Etiqueta LMFDB del carácter (ej: '5.2')  
        max_zeros (int): Número máximo de ceros a recuperar  

    Devuelve:  
        list: Partes imaginarias de ceros (γ)  
    """  
    url = f"https://beta.lmfdb.org/api/lfunctions/degree1/?label={chi_label}&_fields=zeros"  
    try:  
        response = requests.get(url, timeout=10)  
        data = json.loads(response.text)  
        if "zeros" not in data["data"][0]:  
            raise ValueError("Formato LMFDB inválido")  
        return sorted([float(z) for z in data["data"][0]["zeros"] if z != '0'])[:max_zeros]  
    except Exception as e:  
        raise RuntimeError(f"Error LMFDB API: {str(e)}")  

def construir_Hq(q, chi, p_max=100):  
    # ... (igual que en la versión anterior) ...  

def validar_ceros_L(q, chi, tol=1e-5):  
    r"""  
    Versión generalizada con LMFDB.  
    """  
    H = construir_Hq(q, chi)  
    autovalores = sorted([x.real() for x in H.eigenvalues()], key=abs)  

    # Obtener etiqueta LMFDB del carácter  
    n = chi.number()  
    chi_label = f"{q}.{n}"  # Convención LMFDB: q.numero_del_carácter  

    gamma_teoricos = obtener_ceros_LMFDB(q, chi_label)  
    gamma_calculados = [2*pi*q*λ for λ in autovalores[:len(gamma_teoricos)]]  

    datos = {  
        "γ (teórico)": gamma_teoricos,  
        "γ (calculado)": gamma_calculados,  
        "error": [abs(a - b) for a,b in zip(gamma_teoricos, gamma_calculados)]  
    }  
    df = pd.DataFrame(datos)  
    assert all(df["error"] < tol), f"Errores > {tol}: \n{df}"  
    return df  

# --------------------------
# Ejemplo ejecutable (genérico)  
# --------------------------  
if __name__ == "__main__":  
    q = 5  # Cambiar por cualquier módulo  
    G = DirichletGroup(q)  
    chi = G[1]  
    
    print(f"=== Validación HRG para q={q}, χ={chi} ===")  
    try:  
        df = validar_ceros_L(q, chi)  
        print("\nValidación exitosa contra LMFDB:")  
        print(df.to_markdown(index=False))  
    except Exception as e:  
        print(f"\nError: {str(e)}")  