# HRG_Operador.sage
# Versión 2.1 - Código completo con validación LMFDB
# Licencia MIT: https://opensource.org/licenses/MIT

from sage.all import *
import requests
import json
from pandas import DataFrame

def obtener_ceros_LMFDB(q, chi_label, max_zeros=10):
    r"""
    Obtiene ceros no triviales de L(s,χ) desde LMFDB API.
    
    Parámetros:
        q (int): Módulo del carácter
        chi_label (str): Etiqueta LMFDB (ej: '5.2' para q=5, número de carácter 2)
        max_zeros (int): Máximo número de ceros a recuperar
        
    Devuelve:
        list: Partes imaginarias de ceros (γ_n)
    """
    url = f"https://beta.lmfdb.org/api/lfunctions/degree1/?label={chi_label}&_fields=zeros"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()
        zeros = data["data"][0]["zeros"]
        return sorted([float(z) for z in zeros if z != '0'])[:max_zeros]
    except Exception as e:
        raise RuntimeError(f"Error al acceder a LMFDB: {str(e)}")

def construir_Hq(q, chi, p_max=100):
    r"""
    Construye el operador autoadjunto Hq en el espacio modular.
    
    Parámetros:
        q (int): Módulo ≥1
        chi (DirichletCharacter): Carácter primitivo
        p_max (int): Cota superior para primos
        
    Devuelve:
        Matrix: Matriz compleja de dimensión ϕ(q) × ϕ(q)
    """
    if q < 1:
        raise ValueError("El módulo q debe ser ≥1")
    if not chi.is_primitive():
        raise ValueError("Solo caracteres primitivos están implementados")

    G = chi.parent()
    unidades = [k for k in range(1, q) if gcd(k, q) == 1]
    n = len(unidades)
    H = matrix(CDF, n, n, 0.0)
    
    for p in primes(2, p_max):
        if gcd(p, q) != 1:
            continue
            
        try:
            kp = inverse_mod(p, q)
        except:
            continue  # p no es invertible módulo q (improbable si gcd=1)

        for i, m in enumerate(unidades):
            for j, k in enumerate(unidades):
                term = chi(m) * exp(2*pi*CC.0 * m * (k + kp) / p)
                term /= 2 * sqrt(p) * q  # Normalización crítica
                H[i, j] += term

    # Garantizar hermiticidad
    H += H.conjugate().transpose()
    return H

def validar_HRG(q, chi, tol=1e-5):
    r"""
    Valida la correspondencia entre autovalores y ceros de L(s,χ).
    
    Parámetros:
        q (int): Módulo
        chi (DirichletCharacter): Carácter
        tol (float): Tolerancia para errores
        
    Devuelve:
        DataFrame: Resultados comparativos
    """
    print(f"\nConstruyendo Hq para q={q}, χ={chi}...")
    H = construir_Hq(q, chi)
    
    print("Calculando autovalores...")
    autovalores = sorted([x.real() for x in H.eigenvalues() if abs(x.imag()) < 1e-10], key=abs)
    
    print("Obteniendo ceros teóricos de LMFDB...")
    chi_label = f"{q}.{chi.number()}"  # Formato LMFDB
    gamma_teoricos = obtener_ceros_LMFDB(q, chi_label)
    
    gamma_calculados = [2*pi*q*λ for λ in autovalores[:len(gamma_teoricos)]]
    
    datos = {
        "γ (teórico)": gamma_teoricos,
        "γ (calculado)": gamma_calculados,
        "error": [abs(a - b) for a,b in zip(gamma_teoricos, gamma_calculados)]
    }
    df = DataFrame(datos)
    
    if not all(df["error"] < tol):
        print(f"¡Advertencia! Errores superan la tolerancia:\n{df}")
    else:
        print("✓ Validación exitosa dentro de la tolerancia")
    
    return df

# --------------------------
# Ejemplo ejecutable
# --------------------------
if __name__ == "__main__":
    print("=== Demostración de la Correspondencia HRG ===")
    
    # Configuración (cambiar según necesidad)
    q_ejemplo = 5
    G = DirichletGroup(q_ejemplo)
    chi_ejemplo = G[1]  # Segundo carácter (no principal para q=5)
    
    try:
        # Validación completa
        resultados = validar_HRG(q_ejemplo, chi_ejemplo)
        
        # Mostrar resultados en formato tabular
        print("\nResultados finales:")
        print(resultados.to_markdown(index=False))
        
    except Exception as e:
        print(f"\n❌ Error crítico: {str(e)}")
        print("Posibles soluciones:")
        print("1. Verifica tu conexión a internet para acceder a LMFDB")
        print("2. Asegúrate de que SageMath esté correctamente instalado")
        print("3. Comprueba que q y χ son válidos (χ debe ser primitivo)")
