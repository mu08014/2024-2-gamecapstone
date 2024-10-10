Shader "Custom/Kajiya"
{
    Properties
    {
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularTerm("SpecularTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" "Queue"="Geometry" }
        
        Pass
        {
            CGPROGRAM
            #pragma target 3.0
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_shadowcaster
            #include "UnityCG.cginc"
            
            float _DiffuseTerm;
            float _SpecularTerm;
            float _SpecularPower;
            fixed4 _Color;
            sampler2D _MainTex;
            
            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
                float3 normal : TEXCOORD2;
                float2 uv : TEXCOORD3;
            };
            
            struct v2f
            {
                float3 color : COLOR;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float3 normal : NORMAL;

                float4 lightDir : TEXCOORD1;
                float4 cameraDir : TEXCOORD2;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.normal = v.normal;
                o.tangent = normalize(v.tangent);
                
                o.lightDir = mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex);

                o.color = tex2Dlod(_MainTex, float4(v.uv, 0, 0));
                
                return o;
            }
            
            float4 frag (v2f i) : SV_Target
            {
                float diffuse_val = saturate(dot(i.lightDir, i.normal));
                float3 diffuse = i.color;

                float3 halfVec = normalize(i.lightDir + i.cameraDir);
                float3 specular = pow( 
                    sqrt(1 - dot(i.tangent, halfVec) * dot(i.tangent, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                  
                float3 embient = i.color * 0.2f;

                float3 color = (diffuse + specular) * diffuse_val  + embient;
                return float4(color, 1);
            }
            ENDCG
        }

        Pass
        {
            Tags { "LightMode" = "ShadowCaster" }

            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_shadowcaster
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"

            struct v2f
            {
                V2F_SHADOW_CASTER;
            };

            v2f vert(appdata_base v)
            {
                v2f o;
                TRANSFER_SHADOW_CASTER(o)
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                return 0;
            }
            ENDCG
        }

        UsePass "VertexLit/SHADOWCASTER"
    }
    FallBack "Diffuse"
}
